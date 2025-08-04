library(tidyverse)
library(terra)
library(sf)
library(lidR)
library(data.table)
library(minpack.lm)
library(gstat)

# Read in the data that we will need
gediMatched <- st_read("/share/tcsi/lagoodal/GEDI/Storage/gediMatched.shp")
fortypcdBWCA_UTM <- rast("/share/tcsi/lagoodal/GEDI/Storage/fortypcdBWCA.tif") %>% project(., "epsg:26915", method = "near", res = c(30,30))
treeMN <- fread("/share/tcsi/lagoodal/GEDI/Storage/MN_TREE.csv")
cbdRast <- rast("/share/tcsi/lagoodal/GEDI/Storage/cbdRast.tif")
canopyCoverRast <- rast("/share/tcsi/lagoodal/GEDI/Storage/canopyCoverRast.tif")
source("/share/tcsi/lagoodal/GEDI/Storage/lidarMetrics.R")

# Read in and set analysis parameters of the lidar point clouds
LAS <- readLAScatalog("/share/tcsi/lagoodal/GEDI/BWCA_Tiles/", select = "xyzc")
opt_chunk_size(LAS) <- 100
opt_chunk_buffer(LAS) <- 15
LAS@data$code <- str_extract(LAS@data$filename, "[A-Z]{3}[0-9]{6}")

# Create a dataframe of the rows that match the lidar point cloud locations
gediDf <- gediMatched %>% filter(code %in% unique(LAS@data$code))
gediCodes <- unique(gediDf$code)
matchedCodes <- LAS@data$code[LAS@data$code %in% gediCodes]
LAS <- LAS[LAS@data$code %in% gediCodes, ]

# Create folders for metricsDir and alsDir
metricsDir <- "/share/tcsi/lagoodal/GEDI/Tree-Level_Metrics"
alsDir <- "/share/tcsi/lagoodal/GEDI/LAZ"
for (dir in c(metricsDir, alsDir)) if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

# Create the metrics and LAS files.
# TODO parallelize this for the HPC, sequential for now
for (i in 1:nrow(gediDf)) {
  process_tlm2(i, gedi = gediDf, .ctg = LAS, alsdir = alsDir, mdir = metricsDir)
  }  
# Read in the metrics that have ben created by the for loop above
filesMetrics <- list.files(metricsDir, full.names = TRUE)
filesALS <- list.files(alsDir, full.names = TRUE)
idsMetrics <- parse_number(basename(filesMetrics))
idsALS <- parse_number(basename(filesALS))
identical(idsALS, idsMetrics)

# Function to load spatial df, add uniqueID column and rename Zmax to Z
read_with_id <- function(path) {
  id = str_extract(basename(path), "(?<=uID_)[A-Za-z0-9_]+(?=\\.gpkg$)")
  m = st_read(path, quiet = TRUE)
  m$uniqueID = id
  return(select(m, uniqueID, treeID, Z = Zmax, everything()))
}

# Get metrics and then the forest types beneath them
metrics <- map_dfr(filesMetrics, read_with_id) %>% `st_crs<-`("epsg26915")
fortypGEDI <- extract(fortypcdBWCA_UTM, metrics)
metrics <- cbind(metrics, fortypGEDI) %>% rename("Species" = "Label")

# Here I handle the missing species labels
dominantSpp <- metrics %>%
  st_drop_geometry() %>%
  goup_by(uniqueID) %>%
  count(Species) %>%
  summarise(Ndom = sum(n),
            domSpec = Species[which.max(n)],
            percDom = n[which.max(n)] / Ndom,
            Nna = sum(is.na(Species)),
            percNA = Nna/Ndom, .groups = "drop")

metricsDf <- inner_join(metrics, select(dominantSpp, uniqueID, domSpec), by = "uniqueID") %>%
  mutate(Species = case_when(is.na(Species) ~ domSpec, TRUE ~ Species)) %>%
  filter(sum(is.na(Species) == 0, .by = uniquqID))

#----PREDICT DBH----#
# Gotta build the models for predicting DBH. I will use the height (Z) of the trees to
# create models that express a height-diameter relationship
spcds <- c(105, 125, 129, 12, 94, 97, 95, 71, 241, 129, 100, 802, 833, 823, 972, 316, 543, 920, 461, 317, 318, 951, 746, 375, 741, 761)
treeMN <- treeMN %>% filter(SPCD %in% spcds, INVYR == 2023) %>%
  mutate(DIA_CM = DIA * 2.54, HT_M = HT * 0.3048) %>%
  mutate(Species = case_when(SPCD == 12 ~ "Balsam fir", SPCD == 71 ~ "Tamarack", SPCD == 94 ~ "White spruce", SPCD == 95 ~ "Black spruce",
                             SPCD == 97 ~ "Red spruce", SPCD == 105 ~ "Jack pine", SPCD == 125 ~ "Red pine", SPCD == 129 ~ "Eastern white pine",
                             SPCD == 241 ~ "Northern white-cedar", SPCD == 316 ~ "Red maple", SPCD == 317 ~ "Silver maple", SPCD == 318 ~ "Sugar maple",
                             SPCD == 375 ~ "Paper birch", SPCD == 543 ~ "Black ash", SPCD == 741 ~ "Balsam poplar", SPCD == 746 ~ "Quaking aspen",
                             SPCD == 802 ~ "White oak", SPCD == 823 ~ "Bur oak", SPCD == 833 ~ "Northern red oak", SPCD == 951 ~ "American basswood",
                             TRUE ~ "American elm"))

# Create a function to fit Chapman-Richards growth model. I couldn't find any good height-diameter
# relationships in the literature so I decided to estimate them using the nlsLM function.

fitChapman  <- function(df) {
  tryCatch({
    nlsLM(DIA_CM ~ a * (1 - exp(-b * HT_M))^c,
          data = df,
          start = list(a = 30, b = 0.02, c = 1),
          control = nls.lm.control(maxiter = 1000))
  }, error = function(x) {
    # Fallback to log-log linear just incase there is no convergence
    lm(log(DIA_CM) ~ log(HT_M), data = df)
  })
}

fittedSpecies <- treeMN %>% group_by(Species) %>% group_map(~ fitChapman(.x), .keep = TRUE)
species_names <- treeMN %>% group_by(Species) %>% group_keys() %>% pull(Species)

# Because I am trying to assign species to height-diameter curves to forest type code there is
# sometimes a mismatch. So to circumvent this I will assign species as appropriately as possible
# and for those that have more than one (e.g. Eastern white pine / northern red oak / white ash),
# I will assign one species randomly
modelTable <- tibble(
  species = species_names,
  fitModel = fittedSpecies
)

uniqueMetricSpp <- unique(metricsDf$Species)
modelSpp <- unique(modelTable$species)

set.seed(33)
speciesMap <- tibble(
  metricsSpp = uniqueMetricSpp,
  matchedModelSpp = map_chr(unique(metricsDf$Species), function(ms) {
    components <- str_split(tolower(ms), " / |,| and | & ") %>% unlist() %>% str_trim()
    matches <- modelSpp[map_lgl(tolower(modelSpp), function(mod) {
      any(str_detect(mod, fixed(components)))
    })]
    if (length(matches) == 0) {
      NA_character_
    } else {
      sample(matches, 1)
    }
  })
)

# Join with metricsDf and then run the calculation for DBH
metricsDf <- metricsDf %>% left_join(speciesMap, by = c("Species" = "metricsSpp"))
for (i in 1:nrow(modelTable)) {
  S <- modelTable$species[i]
  matchRows <- metricsDf$matchedModelSpp == S
  if (any(matchRows)) {
    metricsDf$DBH[matchRows] <- predict(modelTable$fitModel[[i]], metricsDf[matchRows, ] %>% select(Z)) %>%
      round()
  }
}

#----CROWN CLASS----#
# Function to find crown class
crcl <- function(Z) {
  if (length(Z) == 1) {
    return(c(cl = as.factor("D")))
  } else {
    m <- mean(Z, na.rm = TRUE)
    s <- sd(Z, na.rm = TRUE)
    up <- m + s
    lo <- m - s
    cut(Z, breaks = c(0, lo, up, 100), labels = c("S", "C", "D"))
  }
}
metricsDf$CrCl <- group_by(metricsDf, uniqueID) %>% group_map(~ crcl(.$Z)) %>% unlist()

# Check for any trees that have an unrealistic. If the CBH is too high up, or if
# it is less than 2m. Any tree with an unrealistic CBH is given the default value of Z/2
metricsDf <- metricsDf %>% mutate(CBH = ifelse(CBH > Z * 0.9 | CBH < 2, Z / 2, CBH))

#----CANOPY FUEL LOAD----#
# First we need to get the canopy bulk density values into the metricsDf dataframe
# Get the CBD values located beneath my gedi footprints and add to metricsDf
cbdValues <- extract(cbdRast, metricsDf, ID = FALSE)
metricsDf$CBD <- cbdValues[,1]
metricsDf <- metricsDf %>%
  mutate(CBD = as.character(CBD)) %>%
  mutate(CBD = case_when(
    CBD = "CBD = 1 kg/m^3 X 100" ~ 1 / 100,
    CBD = "CBD = 3 kg/m^3 X 100" ~ 3 / 100,
    CBD = "CBD = 4 kg/m^3 X 100" ~ 4 / 100,
    CBD = "CBD = 5 kg/m^3 X 100" ~ 5 / 100,
    CBD = "CBD = 6 kg/m^3 X 100" ~ 6 / 100,
    CBD = "CBD = 8 kg/m^3 X 100" ~ 8 / 100,
    CBD = "CBD = 9 kg/m^3 X 100" ~ 9 / 100,
    CBD = "CBD = 11 kg/m^3 X 100" ~ 11 / 100,
    CBD = "CBD = 12 kg/m^3 X 100" ~ 12 / 100,
    CBD = "CBD = 16 kg/m^3 X 100" ~ 16 / 100,
    CBD = "CBD = 17 kg/m^3 X 100" ~ 17 / 100,
    CBD = "CBD = 22 kg/m^3 X 100" ~ 22 / 100,
    CBD = "CBD = 24 kg/m^3 X 100" ~ 24 / 100,
    CBD = "CBD = 30 kg/m^3 X 100" ~ 30 / 100,
    CBD = "CBD = 31 kg/m^3 X 100" ~ 31 / 100,
    CBD = "CBD = 34 kg/m^3 X 100" ~ 34 / 100,
    CBD = "CBD = 35 kg/m^3 X 100" ~ 35 / 100,
    TRUE ~ NA_real_
  ))

# Now we can calculate the canopy fuel load by multiplying CBD and Canopy Volume/m2
# CFL units are kb/m2
metricsDf <- metricsDf %>% mutate(CFL = CBD * volume_per_m2)
metricsDf <- st_as_sf(metricsDf)

#----CANOPY COVER----#
ccValues <- extract(canopyCoverRast, metricsDf, ID = FALSE)
metricsDf$CCover <- ccValues[,1]


write_csv(metricsDf, "/share/tcsi/lagoodal/GEDI/Outputs/metricsDf.csv")











































