# prepare rasters for An stephensi modelling

library(terra)
library(tidyverse)

# aggregate all rasters to the same resolution as the nichemapr outputs

# load mechanistically-modelled relative abundance of An stephensi adults per larval
# habitat, aggregate to an annual average abundance, and renormalise
climatic_rel_abund_monthly <- rast("output/rasters/derived/An_stephensi_mechanistic_abundance.tif")
climatic_rel_abund <- mean(climatic_rel_abund_monthly)

# load template raster
covmask <- rast("output/rasters/covariates/covmask.grd")

# align all the covariates to the same resolution and extent as the climatic
# relative abundance layer
agg_ratio <- round(mean(res(climatic_rel_abund) / res(covmask)))
covmask <- aggregate(covmask, agg_ratio)

# resample the climatic relative abundance to make sure it's aligned
climatic_rel_abund <- resample(climatic_rel_abund,
                               covmask,
                               method = "bilinear")

# renormalise this to stretch from 0-1
min <- global(climatic_rel_abund, fun = "min", na.rm = TRUE)[1,1]
climatic_rel_abund <- climatic_rel_abund - min
max <- global(climatic_rel_abund, fun = "max", na.rm = TRUE)[1,1]
climatic_rel_abund <- climatic_rel_abund / max

# mask covmask by this to only use that extent
covmask <- mask(covmask, climatic_rel_abund)

# load all covariates
accessibility <- rast("output/rasters/covariates/accessibility.grd")
built_height <- rast("output/rasters/covariates/built_height.grd")
built_surface <- rast("output/rasters/covariates/built_surface.grd")
built_volume <- rast("output/rasters/covariates/built_volume.grd")
built_c <- rast("output/rasters/covariates/built_c.grd")
smod <- rast("output/rasters/covariates/smod.grd")
landcover <- rast("output/rasters/covariates/landcover.grd")[["landcover_2020"]]
lst_day <- rast("output/rasters/covariates/lst_day.grd")[["lst_day_2020"]]
lst_night <- rast("output/rasters/covariates/lst_night.grd")[["lst_night_2020"]]
nighttimelights <- rast("output/rasters/covariates/nighttimelights.grd")[["nighttimelights_2020"]]
rainfall <- rast("output/rasters/covariates/rainfall.grd")[["rainfall_2020"]]
tcb <- rast("output/rasters/covariates/tcb.grd")[["tcb_2020"]]
tcw <- rast("output/rasters/covariates/tcw.grd")[["tcw_2020"]]
pop <- rast("output/rasters/covariates/pop.grd")[["pop_2020"]]

# set levels for factor rasters
landcover_lookup <- tribble(
  ~id, ~cover,
  00, "Unclassified",
  01, "Evergreen_Needleleaf_Forest",
  02, "Evergreen_Broadleaf_Forest",
  03, "Deciduous_Needleleaf_Forest",
  04, "Deciduous_Broadleaf_Forest",
  05, "Mixed_Forest",
  06, "Closed_Shrublands",
  07, "Open_Shrublands",
  08, "Woody_Savannas",
  09, "Savannas",
  10, "Grasslands",
  11, "Permanent_Wetlands",
  12, "Croplands",
  13, "Urban_And_Built_Up",
  14, "Cropland_Natural_Vegetation_Mosaic",
  15, "Snow_And_Ice",
  16, "Barren_Or_Sparsely_Populated",
  17, "Water"
)
levels(landcover) <- landcover_lookup

smod_lookup <- tribble(
  ~id, ~type,
  30, "URBAN CENTRE",
  23, "DENSE URBAN CLUSTER",
  22, "SEMI-DENSE URBAN CLUSTER",
  21, "SUBURBAN OR PERI-URBAN",
  13, "RURAL CLUSTER",
  12, "LOW DENSITY RURAL",
  11, "VERY LOW DENSITY RURAL",
  10, "WATER"
)
levels(smod) <- smod_lookup

# select, transform, and scale some covariates of larval habitat abundance

larval_covs_continuous <- c(accessibility,
                            log1p(built_volume),
                            log1p(nighttimelights),
                            tcb,
                            tcw,
                            log1p(pop))
larval_covs_continuous <- aggregate(larval_covs_continuous,
                                    agg_ratio,
                                    fun = "mean")
larval_covs_continuous <- scale(larval_covs_continuous)

larval_covs_discrete <- c(smod,
                          landcover)
larval_covs_discrete <- aggregate(larval_covs_discrete,
                                  agg_ratio,
                                  fun = "modal")

larval_covs <- c(larval_covs_continuous,
                 larval_covs_discrete)

# make a new mask based on all of these
new_mask <- app(larval_covs, function(x) !any(is.na(x)))
new_mask[new_mask == 0] <- NA

# tidy up names
names(larval_covs) <- gsub("_2020", "", names(larval_covs))
# and mask to the study area
larval_covs <- mask(larval_covs, new_mask)

# define the populated areas (those to use in the spread model)
populated <- smod %in% c("URBAN CENTRE",
                         "DENSE URBAN CLUSTER",
                         "SEMI-DENSE URBAN CLUSTER",
                         "SUBURBAN OR PERI-URBAN",
                         "RURAL CLUSTER",
                         "LOW DENSITY RURAL")
populated <- aggregate(populated,
                       agg_ratio,
                       fun = "max")
populated <- mask(populated, new_mask)

terra::writeRaster(
  x = new_mask,
  filename = "output/rasters/derived/mask.tif",
  overwrite = TRUE
)

terra::writeRaster(
  x = larval_covs,
  filename = "output/rasters/derived/larval_habitat_covariates.tif",
  overwrite = TRUE
)

terra::writeRaster(
  x = populated,
  filename = "output/rasters/derived/populated_areas.tif",
  overwrite = TRUE
)

terra::writeRaster(
  x = climatic_rel_abund,
  filename = "output/rasters/derived/climatic_rel_abund.tif",
  overwrite = TRUE
)


