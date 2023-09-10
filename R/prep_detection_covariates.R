# make rasters of detection covariates
library(terra)
library(sf)
library(tidyverse)

# raster of populated areas
populated <- rast("~/Dropbox/github/anopheles_stephensi_expansion/output/rasters/derived/populated_areas.tif")

# An. stephensi detections
detections <- readRDS("~/Dropbox/github/anopheles_stephensi_expansion/output/tabular/first_detection.RDS")

# GBIF background point detections
bg <- readRDS("~/Dropbox/github/anopheles_stephensi_expansion/output/tabular/bg_moz_many_20230904.RDS")

# do some reprojection, but stash the originals
mask_old <- populated * 0
old_crs <- "EPSG:4326"
new_crs <- "ESRI:102022"
populated <- terra::project(populated, new_crs)
mask <- populated * 0

bg_observations <- bg %>%
  st_as_sf(
    crs = old_crs,
    coords = c("lon", "lat")
  ) %>%
  st_transform(new_crs) %>%
  mutate(
    populated = terra::extract(populated, .)[, 2],
    valid = !is.na(populated)
  ) %>%
  filter(
    valid & populated
  ) %>%
  sample_n(
    min(100000, n())
  ) %>%
  st_coordinates() %>%
  as.matrix()

plot(populated * 0)
points(bg_observations, cex = 0.5)


bg_count <- terra::rasterize(bg_observations,
                             populated,
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
    populated = terra::extract(populated, .)[, 2],
    valid = !is.na(populated)
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
                     background = 0) %>%
    terra::mask(mask)
  
  subset_density <- focal(subset_count,
                          w = window,
                          fun = "sum",
                          na.policy = "omit",
                          na.rm = TRUE)
  names(subset_density) <- year
  as_detection_density[[match(year, years)]] <- subset_density
}

# project back ot wgs84 and save
as_detection_density_wgs84 <- as_detection_density %>%
  project(old_crs) %>%
  resample(mask_old) %>%
  mask(mask_old)

terra::writeRaster(
  x = as_detection_density_wgs84,
  filename = "output/rasters/derived/an_stephensi_detection_density.tif",
  overwrite = TRUE
)

