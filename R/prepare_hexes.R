library(dplyr)
library(terra)
library(sf)
library(h3) #remotes::install_github("crazycapivara/h3-r")

sapply(
  list.files("R/functions", full.names = TRUE),
  source
)

# get region of interest
source("R/bounding_box.R")

# put hexagons around them
region_hex <- get_h3_from_sf(
  sf_object = st_as_sf(region_shape_buffer))


# filter to only hexes populated areas 
#populated <- rast("output/rasters/derived/populated_areas.tif")
populated <- rast("output/rasters/covariates/populated.tif")


hexes <- region_hex %>%
  sf::st_transform(crs(populated)) %>%
  mutate(
    pop = terra::extract(populated, ., fun = "sum", na.rm = TRUE)[, 2]
  ) %>%
  filter(
    pop > 0
  ) %>%
  mutate(
    id = match(h3_index, unique(h3_index))
  ) %>%
  select(
    -pop
  )

# mask by populated area to make the raster need for zonal stats calculations
hex_raster <- rasterize(hexes,
                        populated,
                        field = "id",
                        fun = min)

hex_raster <- hex_raster * populated

writeRaster(hex_raster, "output/rasters/derived/hex_raster_lookup.tif")
saveRDS(hexes, "output/hexes.RDS")

# compare hex raster with populated raster and region buffer
plot(populated, col = inferno(5)[3:4])
plot(hex_raster, col = viridis::mako(123), add = TRUE)
lines(region_shape_buffer, lwd = 3, col = "darkgoldenrod1")


## alternative hex mapping
region_extent <- ext(populated) %>%
  vect %>%
  st_as_sf()

region_hex_alt <- get_h3_from_sf(region_extent)

hexes_alt <- region_hex_alt %>%
  sf::st_transform(crs(populated)) %>%
  mutate(
    pop = terra::extract(populated, ., fun = "sum", na.rm = TRUE)[, 2]
  ) %>%
  filter(
    pop > 0
  ) %>%
  mutate(
    id = match(h3_index, unique(h3_index))
  ) %>%
  select(
    -pop
  )

# mask by populated area to make the raster need for zonal stats calculations
hex_raster_alt <- rasterize(hexes_alt,
                        populated,
                        field = "id",
                        fun = min)

hex_raster_alt <- hex_raster_alt * populated


# compare hex raster with populated raster and region buffer
plot(populated, col = inferno(5)[3:4])
plot(hex_raster_alt, col = viridis::plasma(123), add = TRUE)
lines(region_shape_buffer, lwd = 3, col = "darkgoldenrod1")

# compare hexes and alternative hexes
plot(populated, col = inferno(5)[3:4])
plot(hex_raster_alt, col = viridis::plasma(123), add = TRUE)
plot(hex_raster, col = viridis::mako(123), add = TRUE)
