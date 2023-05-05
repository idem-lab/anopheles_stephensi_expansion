library(terra)
library(dplyr)

##### accessibility
# Accessibility to cities for a nonimal year 2015.
# "Cities" are defined as contiguous areas with 1,500
# or more inhabitants per square kilometre or a majority
# of built-up land cover types coincident with a 
# population centre of at least 50,000 inhabitants. Pixel
# values show estimated fasted land-based travel time to
# the nearest city in minutes. Produced by Dr Dan Weiss
# (https://doi.org/10.1038/nature25181).

accessibility <- rast("data/MAP_covariates/Accessibility/accessibility_to_cities_2015_v1.0.tif")
names(accessibility) <- "accessibility"

accessibility

plot(accessibility)


#### EVI 
# EVI is derived from the 8-daily global 1km MODIS
# v6 MCD43D62, MCD43D63 and MCD43D64 products. 
# This is then gapfilled using an algorithm 
# developed by Dr Dan Weiss and implemented 
# globally by Dr Harry Gibson 
# (https://doi.org/10.1016/j.isprsjprs.2014.10.001). 
# The gapfilled outputs are aggregated temporally
# to the annual level using a mean. 


evi <- list.files(
  path = "data/MAP_covariates/evi/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

evi_names <- names(evi) %>%
  sub(
    pattern = ".*v6\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("evi_", .)

names(evi) <- evi_names

evi

plot(c(evi$evi_2015, evi$evi_2021))

#### GHS Built
# Built-up area grids  derived from the Global Land 
# Survey (GLS) Landsat image collections (GLS1975, 
# GLS1990, GLS2000, and ad-hoc Landsat 8 collection 
# 2013/2014) (https://ghsl.jrc.ec.europa.eu/ghs_bu2019.php). 
# Pixel values represent built-up area density from 0-100,
# aggregated from the 30m Landsat data.

ghs_built <- rast("data/MAP_covariates/GHSL_old/GHS_BUILT_R18A_v2.2014.Annual.Data.1km.Data.tif")

names(ghs_built) <- "ghs_built"

ghs_built

plot(ghs_built)

#### GHS SMOD
# Settlement grids delineating and classifying settlement
# typologies via a logic of population size, population
# and built-up area densities 
# (https://ghsl.jrc.ec.europa.eu/ghs_smod2019.php). 
# The pixel classification criteria are available in the
# supporting data package PDF.

ghs_smod <- rast("data/MAP_covariates/GHSL_old/GHS_SMOD_R19A_v2.2015.Annual.Data.1km.Data.tif")

names(ghs_smod) <- "ghs_smod"

ghs_smod

plot(ghs_smod)

### Landcover
# Landcover classification data derived from MODIS v6
# MCD12Q1, using the IGBP classification. Annual
# majority rasters (class number covering the majority
# of each pixel) are available, derived and aggregated
# from the 500m original datasets.
# 
# IGBP Landcover Classes:
#   00 Unclassified
# 01 Evergreen_Needleleaf_Forest
# 02 Evergreen_Broadleaf_Forest
# 03 Deciduous_Needleleaf_Forest
# 04 Deciduous_Broadleaf_Forest
# 05 Mixed_Forest
# 06 Closed_Shrublands
# 07 Open_Shrublands
# 08 Woody_Savannas
# 09 Savannas
# 10 Grasslands
# 11 Permanent_Wetlands
# 12 Croplands
# 13 Urban_And_Built_Up
# 14 Cropland_Natural_Vegetation_Mosaic
# 15 Snow_And_Ice
# 16 Barren_Or_Sparsely_Populated
# 17 Water

landcover <- list.files(
  path = "data/MAP_covariates/Landcover/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

landcover_names <- names(landcover) %>%
  sub(
    pattern = ".*Landcover\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("landcover_", .)

names(landcover) <- landcover_names

landcover

plot(c(landcover$landcover_2015, landcover$landcover_2020))

#### LST Day
# Land surface temperature
# LST_Day is derived from the 8-daily global 1km
# MODIS MOD11A2 v6 products. This is then
# gapfilled using an algorithm developed by Dr 
# Dan Weiss and implemented globally by Dr Harry
# Gibson 
# (https://doi.org/10.1016/j.isprsjprs.2014.10.001). 
# The gapfilled outputs are aggregated
# temporally to the annual level using a mean.

lst_day <- list.files(
  path = "data/MAP_covariates/LST_Day/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

lst_day_names <- names(lst_day) %>%
  sub(
    pattern = ".*v6\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("lst_day_", .)

names(lst_day) <- lst_day_names

lst_day

plot(c(lst_day$lst_day_2015, lst_day$lst_day_2021))

#### LST NIGHT
# LST_NIGHT is derived from the 8-daily global 1km
# MODIS MOD11A2 v6 products. This is then
# gapfilled using an algorithm developed by Dr 
# Dan Weiss and implemented globally by Dr Harry
# Gibson 
# (https://doi.org/10.1016/j.isprsjprs.2014.10.001). 
# The gapfilled outputs are aggregated
# temporally to the annual level using a mean.

lst_night <- list.files(
  path = "data/MAP_covariates/LST_Night/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

lst_night_names <- names(lst_night) %>%
  sub(
    pattern = ".*v6\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("lst_night_", .)

names(lst_night) <- lst_night_names

lst_night

plot(c(lst_night$lst_night_2015, lst_night$lst_night_2021))

#### Night time Lights

nighttimelights <- list.files(
  path = "data/MAP_covariates/Nighttime_Lights/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

nighttimelights_names <- names(nighttimelights) %>%
  sub(
    pattern = ".*Background\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("nighttimelights_", .)

names(nighttimelights) <- nighttimelights_names

nighttimelights <- mask(nighttimelights, ghs_smod)

nighttimelights

plot(c(nighttimelights$nighttimelights_2015, nighttimelights$nighttimelights_2021))

#### rainfall 
# Annual rainfall totals from the CHIRPS dataset
# (https://www.chc.ucsb.edu/data/chirps).
# The 1km version here is a neareast-neighbour
# resample of the lower resolution data
# available from CHIRPS.

rainfall <- list.files(
  path = "data/MAP_covariates/Rainfall/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

rainfall_names <- names(rainfall) %>%
  sub(
    pattern = ".*2-0\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("rainfall_", .)

names(rainfall) <- rainfall_names

rainfall

plot(c(rainfall$rainfall_2015, rainfall$rainfall_2021))

# TCB
# tasselated cap brightness
# TCB is derived from the 8-daily global 1km MODIS
# v6 MCD43D62, MCD43D63 and MCD43D64 products. 
# This is then gapfilled using an algorithm 
# developed by Dr Dan Weiss and implemented 
# globally by Dr Harry Gibson 
# (https://doi.org/10.1016/j.isprsjprs.2014.10.001). 
# The gapfilled outputs are aggregated temporally
# to the annual level using a mean. 


tcb <- list.files(
  path = "data/MAP_covariates/TCB/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

tcb_names <- names(tcb) %>%
  sub(
    pattern = ".*v6\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("tcb_", .)

names(tcb) <- tcb_names

tcb

plot(c(tcb$tcb_2015, tcb$tcb_2021))


# TCW
# tasselated cap wetness
# TCW is derived from the 8-daily global 1km MODIS
# v6 MCD43D62, MCD43D63 and MCD43D64 products. 
# This is then gapfilled using an algorithm 
# developed by Dr Dan Weiss and implemented 
# globally by Dr Harry Gibson 
# (https://doi.org/10.1016/j.isprsjprs.2014.10.001). 
# The gapfilled outputs are aggregated temporally
# to the annual level using a mean. 


tcw <- list.files(
  path = "data/MAP_covariates/TCW/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

tcw_names <- names(tcw) %>%
  sub(
    pattern = ".*v6\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("tcw_", .)

names(tcw) <- tcw_names

tcw

plot(c(tcw$tcw_2015, tcw$tcw_2021))


# Worldpop
# Annual 1km UN-adjusted population counts 
# from WorldPop v3 
# (https://www.worldpop.org/geodata/listing?id=75). 
# This version has been derived by mosaicing the
# country outputs and aligning to MAP's master 
# coastline template (reallocating population from 
# cells falling outside the MAP coastline into the
# nearest land pixel). 

pop <- list.files(
  path = "data/MAP_covariates/WorldPop/",
  full.names = TRUE
) %>%
  sapply(
    FUN = rast
  ) %>%
  rast

pop_names <- names(pop) %>%
  sub(
    pattern = ".*fix\\.",
    replacement = "",
    x = .
  ) %>%
  sub(
    pattern = "\\.Annual.*",
    replacement = "",
    x = .
  ) %>%
  paste0("pop_", .)

names(pop) <- pop_names

pop

plot(c(pop$pop_2015, pop$pop_2020))

# leave out lst, evi, tcb, tcw, rainfall
# pca of others, use top 3 axes, create layers
# 2015 ftw
# 

accessibility
#evi
ghs_built
ghs_smod
landcover
#lst_day
#lst_night
nighttimelights
#rainfall
#tcb
#tcw
pop

covs <- c(
  accessibility,
  ghs_built,
  ghs_smod,
  landcover$landcover_2015,
  nighttimelights$nighttimelights_2015,
  pop$pop_2015
)

covs_continuous <- c(
  accessibility,
  ghs_built,
  #ghs_smod,
  #landcover$landcover_2015,
  nighttimelights$nighttimelights_2015,
  pop$pop_2015
)

cov_dat <- as.matrix(covs_continuous)

cov_dat <- cov_dat[!is.nan(cov_dat[,1]),]

cov_dat <- cov_dat[!is.nan(cov_dat[,4]),] # ~500k values in pop are NaN - explore this

head(cov_dat)
dim(cov_dat)

pca_cov <- prcomp(cov_dat, scale. = TRUE)

summary(pca_cov)
plot(pca_cov)


