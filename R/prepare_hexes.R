library(dplyr)
library(terra)
library(sf)
library(h3) #remotes::install_github("crazycapivara/h3-r")
#library(movement)

sapply(
  list.files("R/functions", full.names = TRUE),
  source
)

# get region of interest
source("R/bounding_box.R")

# put hexagons around them
region_hex <- get_h3_from_sf(
  sf_object = st_as_sf(region_shape_buffer))

plot(region_hex)

# get regions with significant populations 
populated <- rast("output/rasters/covariates/populated.grd")


populated_agg <- aggregate(
  x = populates,
  fact = 5,
  fun = "any"
)

# convert hexes to raster of this population
hexvec <- vect(region_hex) #%>%
  #terra::project("EPSG:4326")

# check if each area is populated
region_hex_pop <- terra::zonal(
  populated_agg,
  hexvec,
  fun = sum,
  na.rm=TRUE
)

#subsets to those with finite population count

populated_hexes <- region_hex[!is.nan(region_hex_pop[,1]),]
plot(populated_hexes)


saveRDS(
  populated_hexes
)
