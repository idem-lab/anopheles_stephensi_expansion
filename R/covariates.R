library(terra)
library(dplyr)

sapply(
  list.files("R/functions/", full.names = TRUE),
  source
)

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

# plot(c(pop$pop_2015, pop$pop_2020))

### GHS_BUILT_H
# Average of the Gross Building Height (AGBH) and Average
# of the Net Building Height (ANBH) for 2018 from GHSL 
# (https://ghsl.jrc.ec.europa.eu/ghs_buH2023.php). Pixel
# values are average height of the built surfaces in 
# meters. The versions here have been aggregated from the
# 100m originals first using a mean in the original 
# mollweide projection, and then reporjected to wgs84 
# using bilinear resampling.

# here using gross built height
built_height <- rast(
  x = "data/MAP_covariates/GHSL_2023/GHS_BUILT_H_AGBH_R23A.2018.Annual.Data.1km.mean.tif"
)

names(built_height) <- "built_height"

built_height

## Make mask from pop layer as this has most NAs - has the lakes and oceans cut out

covmask <- mask(pop$pop_2020, built_height)
names(covmask) <- "mask"
covmask[!is.na(covmask)] <- 0

covmask <- writereadrast(
  covmask,
  "output/rasters/covariates/covmask.grd"
)

covmask
# plot(covmask)

built_height <- terra::mask(built_height, covmask)

built_height <- writereadrast(
  built_height,
  "output/rasters/covariates/built_height.grd"
)

pop <- terra::mask(pop, covmask)

pop <- writereadrast(
  pop,
  "output/rasters/covariates/pop.grd"
)

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

accessibility <- terra::mask(accessibility, covmask)

accessibility <- writereadrast(
  accessibility,
  "output/rasters/covariates/accessibility.grd"
)

accessibility

# plot(accessibility)


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
  path = "data/MAP_covariates/EVI/",
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

evi <- terra::mask(evi, covmask)

evi <- writereadrast(
  evi,
  "output/rasters/covariates/evi.grd"
)

evi

# plot(c(evi$evi_2015, evi$evi_2021))

# GHSL OLD
# 
# #### GHS Built
# # Built-up area grids  derived from the Global Land 
# # Survey (GLS) Landsat image collections (GLS1975, 
# # GLS1990, GLS2000, and ad-hoc Landsat 8 collection 
# # 2013/2014) (https://ghsl.jrc.ec.europa.eu/ghs_bu2019.php). 
# # Pixel values represent built-up area density from 0-100,
# # aggregated from the 30m Landsat data.
# 
# ghs_built <- rast("data/MAP_covariates/GHSL_old/GHS_BUILT_R18A_v2.2014.Annual.Data.1km.Data.tif")
# 
# names(ghs_built) <- "ghs_built"
# 
# ghs_built <- mask(ghs_built, covmask)
# 
# ghs_built
# 
# # plot(ghs_built)
# 
# #### GHS SMOD
# # Settlement grids delineating and classifying settlement
# # typologies via a logic of population size, population
# # and built-up area densities 
# # (https://ghsl.jrc.ec.europa.eu/ghs_smod2019.php). 
# # The pixel classification criteria are available in the
# # supporting data package PDF.
# 
# ghs_smod <- rast("data/MAP_covariates/GHSL_old/GHS_SMOD_R19A_v2.2015.Annual.Data.1km.Data.tif")
# 
# names(ghs_smod) <- "ghs_smod"
# 
# ghs_smod <- mask(ghs_smod, covmask)
# 
# ghs_smod
# 
# # plot(ghs_smod)

## GHSL 2023

### GHS_BUILT_H
# Average of the Gross Building Height (AGBH) and Average
# of the Net Building Height (ANBH) for 2018 from GHSL 
# (https://ghsl.jrc.ec.europa.eu/ghs_buH2023.php). Pixel
# values are average height of the built surfaces in 
# meters. The versions here have been aggregated from the
# 100m originals first using a mean in the original 
# mollweide projection, and then reporjected to wgs84 
# using bilinear resampling.

# here using gross built height
built_height <- rast(
  x = "data/MAP_covariates/GHSL_2023/GHS_BUILT_H_AGBH_R23A.2018.Annual.Data.1km.mean.tif"
) %>%
  terra::mask(
    mask = covmask
  )

names(built_height) <- "built_height"

built_height <- writereadrast(
  built_height,
  "output/rasters/covariates/built_height.grd"
)

built_height

# plot(built_height)

### GHS_BUILT_S
# Built-up surface grid for 2020 from GHSL, for total
# residential and non-residential 
# (https://ghsl.jrc.ec.europa.eu/ghs_buS2023.php). Pixel
# values are built square meters in the grid cell. The
# version here has been reprojected from the 1km
# mollweide dataset to wgs84 using bilinear resampling.
built_surface <- rast(
  x = "data/MAP_covariates/GHSL_2023/GHS_BUILT_S_R23A.2020.Annual.Data.1km.Data.tif"
) %>%
  terra::mask(
    mask = covmask
  )

names(built_surface) <- "built_surface"

built_surface <- writereadrast(
  built_surface,
  "output/rasters/covariates/built_surface.grd"
)

built_surface

# plot(built_surface)

### GHS_BUILT_V
# Built-up volume grid for 2020 from GHSL, for total
# residential and non-residential
# (https://ghsl.jrc.ec.europa.eu/ghs_buV2023.php). 
# Pixel values are built cubic meters in the grid cell. 
# The version here has been reprojected from the 1km
# mollweide dataset to wgs84 using bilinear resampling.
built_volume <- rast(
  x = "data/MAP_covariates/GHSL_2023/GHS_BUILT_V_R23A.2020.Annual.Data.1km.Data.tif"
) %>%
  terra::mask(
    mask = covmask
  )

names(built_volume) <- "built_volume"

built_volume <- writereadrast(
  built_volume,
  "output/rasters/covariates/built_volume.grd"
)

built_volume

# plot(built_volume)


### GHS_SMOD
# Settlement grids delineating and classifying 
# settlement typologies via a logic of population size,
# population and built-up area densities
# (https://ghsl.jrc.ec.europa.eu/ghs_smod2023.php). The
# pixel classification criteria are available in the
# supporting data package PDF. The version here has been
# reprojected from the original 1km mollweide dataset
# using nearest neighbour resampling.
# 30: URBAN CENTRE GRID CELL
# 23: DENSE URBAN CLUSTER GRID CELL
# 22: SEMI-DENSE URBAN CLUSTER GRID CELL
# 21: SUBURBAN OR PERI-URBAN GRID CELL
# 13: RURAL CLUSTER GRID CELL
# 12: LOW DENSITY RURAL GRID CELL
# 11: VERY LOW DENSITY RURAL GRID CELL
# 10: WATER GRID CELL
smod <- rast(
  x = "data/MAP_covariates/GHSL_2023/GHS_SMOD_R23A.2020.Annual.Data.1km.Data.tif"
) %>%
  terra::mask(
    mask = covmask
  )

names(smod) <- "smod"

smod <- writereadrast(
  smod,
  "output/rasters/covariates/smod.grd"
)

smod

# plot(smod)


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

landcover <- terra::mask(landcover, covmask)

landcover <- writereadrast(
  landcover,
  "output/rasters/covariates/landcover.grd"
)

landcover

# plot(c(landcover$landcover_2015, landcover$landcover_2020))

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

lst_night <- mask(lst_day, covmask)

lst_day <- writereadrast(
  lst_day,
  "output/rasters/covariates/lst_day.grd"
)

lst_day

# plot(c(lst_day$lst_day_2015, lst_day$lst_day_2021))

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

lst_night <- mask(lst_night, covmask)

lst_night <- writereadrast(
  lst_night,
  "output/rasters/covariates/lst_night.grd"
)

lst_night

# plot(c(lst_night$lst_night_2015, lst_night$lst_night_2021))

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

nighttimelights <- terra::mask(nighttimelights, covmask)

nighttimelights <- writereadrast(
  nighttimelights,
  "output/rasters/covariates/nighttimelights.grd"
)

nighttimelights

# plot(c(nighttimelights$nighttimelights_2015, nighttimelights$nighttimelights_2021))

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

rainfall <- mask(rainfall, covmask)

rainfall <- writereadrast(
  rainfall,
  "output/rasters/covariates/rainfall.grd"
)

rainfall

# plot(c(rainfall$rainfall_2015, rainfall$rainfall_2021))

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

tcb <- terra::mask(tcb, covmask)

tcb <- writereadrast(
  tcb,
  "output/rasters/covariates/tcb.grd"
)

tcb

# plot(c(tcb$tcb_2015, tcb$tcb_2021))


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

tcw  <- terra::mask(tcw, covmask)

tcw <- writereadrast(
  tcw,
  "output/rasters/covariates/tcw.grd"
)

tcw

# plot(c(tcw$tcw_2015, tcw$tcw_2021))

## ----Unique Human Settlement Boundaries----

# This dataset represents the delineation of unique
# human settlements within the region of interest.
# Using percolation theory on the road network, 
# which are the primary points of interaction and
# human agglomeration, we identify unique human
# settlements. 
# 
# In this case, for the clustering procedure, given
# a graph of the road network, where nodes represent
# intersections and the weight for each edge is the
# length of the street that connects them and a
# certain metric threshold (in this case 500m) we
# produce a network percolation via the following
# steps:
# 1. We select a random node of the graph,
# generating a new cluster and inserting the node
# into the cluster.
# 2. We keep a first-in first-out queue of nodes to
# expand, from which we extract a node to continue
# the process. We add the node in step 1 to this
# queue. Nodes are only added to this queue if they
# are not already included.
# 3. We extract a node from the queue of nodes to
# explore and if a link departing from that node
# (not yet included in the cluster) is smaller than
# the threshold, include the node in the cluster
# and the end node of the link in the queue of
# nodes to explore.
# 4. We repeat step 3 until no further node can be
# expanded (the queue is empty) and if there are
# nodes left in the graph that do not belong to any
# cluster, generate a new cluster by a random
# available node and repeat from step 1.
# 
# Once we have all nodes associated with a unique
# cluster we rasterize the clusters at ~1Km pixel
# size. Where the pixel value is an id of the
# cluster.


# settlement <- rast(
#   x = "data/MAP_covariates/Unique_Human_Settlement_Boundaries/Unique_Settlement_Boundaries_Percolation_1km.tif"
# )
# # the extent of this variable does not match the others
# 
# names(settlement) <- "settlement"
# 
# settlement[] <- ifelse(!is.na(as.array(settlement)), 1, NA) # something wrong
# # plot(settlement)

## PCA of covariates

# leave out lst, evi, tcb, tcw, rainfall
# pca of others, use top 3 axes, create layers
# using 2020 
# 

accessibility
#evi
#ghs_built
#ghs_smod
built_height
built_surface
built_volume
landcover
#lst_day
#lst_night
nighttimelights
#rainfall
#tcb
#tcw
pop

gc()


covs <- c(
  accessibility,
  built_height,
  built_surface,
  built_volume,
  smod,
  landcover$landcover_2020,
  nighttimelights$nighttimelights_2020,
  pop$pop_2020
)

covs <- writereadrast(
  covs,
  "output/rasters/covariates/covs.grd"
)

covs_continuous <- c(
  accessibility,
  built_height,
  built_surface,
  built_volume,
  nighttimelights$nighttimelights_2020,
  pop$pop_2020
)

covs_continuous <- writereadrast(
  covs_continuous,
  "output/rasters/covariates/covs_continuous.grd"
)

covs_continuous <- rast("output/rasters/covariates/covs_continuous.grd")

gc()


prcdat <- as.matrix(covs_continuous)
prcdat <- prcdat[!is.na(prcdat[,1]),]

#prcdat <- covs_continuous[]

gc()

apply(prcdat, 2, function(x){sum(is.na(x))})

samples <- sample(
  x = 1:dim(prcdat)[1],
  size = 1e5,
  replace = FALSE
)

prccovs <- prcomp(prcdat[samples,], scale. = TRUE)

prccovs
summary(prcovs)
# this essenially shows that pc2 is accessibilty, pc3 is nighttime lights, pc4 is population


prbuilt <- prcomp(prcdat[samples, 2:4], scale. = TRUE)

prbuilt
summary(prbuilt)
# 96% coverage by one component. 
# suggest using just built_volume instead of a pc of the three, as volume includes height and surface

pairs(
  prcdat[samples, 2:4]
)

covs <- c(
  accessibility,
  #built_height,
  #built_surface,
  built_volume,
  smod,
  landcover$landcover_2020,
  nighttimelights$nighttimelights_2020,
  pop$pop_2020
)

covs <- writereadrast(
  covs,
  "output/rasters/covariates/covs.grd"
)


# 

covmask <- rast("output/rasters/covariates/covmask.grd")
covs <- rast("output/rasters/covariates/covs.grd")
covs_continuous <- rast("output/rasters/covariates/covs_continuous.grd")

accessibility <- rast("output/rasters/covariates/accessibility.grd")
#evi <- rast("output/rasters/covariates/")
#ghs_built
#ghs_smod
built_height <- rast("output/rasters/covariates/built_height.grd")
built_surface <- rast("output/rasters/covariates/built_surface.grd")
built_volume <- rast("output/rasters/covariates/built_volume.gri")
smod <- rast("output/rasters/covariates/smod.grd")
landcover <- rast("output/rasters/covariates/landcover.grd")
lst_day <- rast("output/rasters/covariates/lst_day.grd")
lst_night <- rast("output/rasters/covariates/lst_night.grd")
nighttimelights <- rast("output/rasters/covariates/nighttimelights.grd")
rainfall <- rast("output/rasters/covariates/rainfall.grd")
tcb <- rast("output/rasters/covariates/tcb.grd")
tcw <- rast("output/rasters/covariates/tcw.grd")
pop <- rast("output/rasters/covariates/pop.grd")
