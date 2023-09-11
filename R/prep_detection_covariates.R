# make rasters of detection covariates
library(terra)
library(sf)
library(tidyverse)

# get combined mask
hex_lookup <- rast("output/rasters/derived/hex_raster_lookup.tif")
mask_old <- rast("output/rasters/derived/mask.tif")
mask_old <- mask(mask_old, hex_lookup)

# An. stephensi detections
detections <- readRDS("output/tabular/first_detection.RDS")

# GBIF background point detections
bg <- readRDS("output/tabular/bg_moz_many_20230904.RDS")

# do some reprojection, but stash the originals
old_crs <- "EPSG:4326"
new_crs <- "ESRI:102022"
mask <- terra::project(mask_old, new_crs)

bg_observations <- bg %>%
  st_as_sf(
    crs = old_crs,
    coords = c("lon", "lat")
  ) %>%
  st_transform(new_crs) %>%
  mutate(
    mask = terra::extract(mask, .)[, 2],
    valid = !is.na(mask)
  ) %>%
  filter(
    valid & mask
  ) %>%
  sample_n(
    min(100000, n())
  ) %>%
  st_coordinates() %>%
  as.matrix()

plot(mask)
points(bg_observations, cex = 0.5)

bg_count <- terra::rasterize(bg_observations,
                             mask,
                             fun = length,
                             background = 0)
bg_count <- mask(bg_count, mask)

# do 2D KDE on this
radius_km <- 250
window <- focalMat(mask,
                   d = radius_km * 1e3,
                   type = "Gauss")

bg_density <- focal(bg_count,
                    w = window,
                    fun = "sum",
                    na.policy = "omit",
                    na.rm = TRUE)


# Now do time-varying An. stephensi detection probabilities

as_detections <- detections %>%
  filter(
    ever_detected == 1
  ) %>%
  st_as_sf(
    crs = old_crs,
    coords = c("x", "y")
  ) %>%
  st_transform(
    new_crs
  ) %>%
  mutate(
    mask = terra::extract(mask, .)[, 2],
    valid = !is.na(mask)
  ) %>%
  filter(
    valid
  ) %>%
  bind_cols(
    st_coordinates(.)
  ) %>%
  select(
    X, Y,
    year_first_detected
  ) %>%
  st_drop_geometry()

years <- seq(min(as_detections$year_first_detected),
             max(as_detections$year_first_detected),
             by = 1)

as_detection_density <- replicate(
  length(years),
  mask,
  simplify = FALSE
) %>%
  do.call(c, .)

for (year in years) {
  
  # for each year, use the data up to the previous year (except in first year in
  # timeseries)
  min_year <- max(year - 1, min(year))
  subset_count <- as_detections %>%
    filter(
      year_first_detected <= min_year
    ) %>%
    select(X, Y) %>%
    as.matrix() %>%
    terra::rasterize(mask,
                     fun = length,
                     background = 0)
  
  subset_density <- focal(subset_count,
                          w = window,
                          fun = "sum",
                          na.policy = "omit",
                          na.rm = TRUE)
  names(subset_density) <- year
  as_detection_density[[match(year, years)]] <- subset_density
}

# project back to wgs84, mask and save
as_detection_density_wgs84 <- as_detection_density %>%
  project(old_crs) %>%
  resample(mask_old) %>%
  mask(mask_old)

terra::writeRaster(
  x = as_detection_density_wgs84,
  filename = "output/rasters/derived/an_stephensi_detection_density.tif",
  overwrite = TRUE
)

