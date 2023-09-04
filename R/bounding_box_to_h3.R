#get matching h3 geometries for bounding box area
library(h3); library(sf); library(terra)
h3_projection <- "ESRI:102022"

source("R/bounding_box.R")


region_shape_buffer_sf <- st_as_sf(region_shape_buffer)

region_hex <- polyfill(region_shape_buffer_sf, res = 1)
region_hex <- unique(region_hex)
region_hex <- h3_to_geo_boundary_sf(region_hex)

region_hex <- st_transform(region_hex,crs = h3_projection)
plot(region_hex)

