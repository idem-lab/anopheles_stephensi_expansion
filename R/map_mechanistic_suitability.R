# Map suitability for An. stephensi and An. gambiae based on Nichemapr
# microclimate model, a cone model of larval habitat, estimated lifehistory
# relationships, and synoptic annual climatic conditions

library(NicheMapR)
library(tidyverse)

source("R/functions/micro_functions.R")

# load the functions giving lifehistory parameters as a function of temperature
lifehistory_functions_stephensi <- get_lifehistory_functions("An. stephensi")
# lifehistory_functions_gambiae <- get_lifehistory_functions("An. gambiae")


# find all the centroids of the climate raster, for the Afro and EMRO regions,
# so we can evaluate this all over the raster

# find and load the global_climate raster downloaded by NicheMapr
gcfolder <- paste(.libPaths()[1], "/gcfolder.rda", sep = "")
load(gcfolder)
nichemapr_climate_raster_mask <- terra::rast(paste0(folder, "/global_climate.nc"))[[1]] * 0 
names(nichemapr_climate_raster_mask) <- "mask"
# crop down to EMRO and AFRO regions

# get a shapefile of the region of interest
library(tidyverse)
library(terra)
library(geodata)
world_country_shapes <- world(level = 0, path = "~/build")
keep <- world_country_shapes$GID_0 %in% region_countries()
region_country_shapes <- world_country_shapes[keep, ]
region_shape <- terra::aggregate(region_country_shapes)

# buffer the region slightly (10km), and use it to create a raster of the
# relevant (non-NA) pixels
region_shape_buffer <- terra::buffer(region_shape, 1e4)

region_raster_mask <- region_shape_buffer %>%
  terra::rasterize(
    nichemapr_climate_raster_mask
  ) %>%
  crop(
    region_shape_buffer
  ) * crop(
    nichemapr_climate_raster_mask,
    region_shape_buffer
  )

# pull out the non-NA cells in this
coords <- xyFromCell(region_raster_mask, cells(region_raster_mask))

# aggregate it to get a lower resolution set of coordinates for debugging -
# original resolution is ~18.5km at the equator
region_raster_mask_agg <- aggregate(region_raster_mask, 15)
coords_agg <- xyFromCell(region_raster_mask_agg, cells(region_raster_mask_agg))
nrow(coords_agg)

# find and skip NAs in the global climate raster at these locations (use raster,
# as terra gives slightly different locations?)
nichemapr_climate_raster <- raster::brick(paste0(folder, "/global_climate.nc"))
vals <- raster::extract(nichemapr_climate_raster, coords_agg)
fine <- apply(is.finite(vals), 1, all)
valid_cells <- which(fine)

# get pixels as a list
coords_list <- coords_agg[valid_cells, ] %>%
  as_tibble() %>%
  split(seq_len(nrow(.)))

# run suitability for all pixels

# for some reason, nichemapr seems to error with any future plan except
# sequential, so go old school and use snowfall!
library(snowfall)
sfInit(parallel = TRUE, cpus = 8)
sfLibrary(NicheMapR)
sfLibrary(tidyverse)
# manually export every required function to the clusters
sfExport(list = list("model_climatic_conditions",
                     "smooth_rainfall",
                     "ensure_positive",
                     "enforce_max",
                     "cone_depth_to_volume",
                     "cone_depth_to_radius",
                     "cone_volume_to_surface",
                     "evaporation_rate_kg_h",
                     "saturation_vapour_pressure_kpa",
                     "air_pressure_kpa",
                     "iterate_cone_volume",
                     "simulate_ephemeral_habitat",
                     "simulate_population",
                     "iterate_state",
                     "create_matrix",
                     "summarise_dynamics"))
time_agg <- system.time(
  results_list <- sfLapply(coords_list,
                           calculate_suitability,
                           lifehistory_functions = lifehistory_functions_stephensi)
)

# # look for the bad ones - could do this faster in parallel with a trycatch
# for (i in seq_along(coords_list)) {
#   print(i)
#   tmp <- model_climatic_conditions(coords_list[[i]])
# }
# 
# model_climatic_conditions(coords_list[[94]])

# # site long 35.3333333333333 lat 32
# 
# dodgy <- which(near(coords_agg[valid_cells, "x"], 35.3333333333333) &
#         near(coords_agg[valid_cells, "y"], 32))
# # dodgy <- 2728
# 
# tmp <- calculate_stephensi_suitability(coords_list[[dodgy]])
# 
# nrow(coords_agg)

# time in minutes - 41 for aggregated
# in parallel, 191 mins @ agg 4
time_agg["elapsed"] / 60

sfStop()

# add on the coordinates
index_list <- lapply(seq_along(coords_list), function(x) tibble(cell_index = x))

results <- mapply(bind_cols, results_list, coords_list, index_list, SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  `rownames<-`(NULL)

head(results)
tail(results)

results_annual <- results %>%
  filter(
    microclimate == "habitat",
    larval_habitat == "permanent"
  ) %>%
  group_by(
    cell_index,
    x,
    y,
    microclimate
  ) %>%
  summarise(
    months_suitable = sum(persistence),
    relative_abundance = mean(relative_abundance),
    .groups = "drop"
  )
# mutate(
#   relative_abundance = case_when(
#     months_suitable == 0 ~ 0,
#     .default = relative_abundance
#   )
# )


months_suitable_agg <- annual_relative_abundance_agg <- region_raster_mask_agg
names(months_suitable_agg) <- "months suitable"
names(annual_relative_abundance_agg) <- "relative abundance (annual)"

cell_index <- cellFromXY(region_raster_mask_agg, coords_agg[valid_cells, ])
months_suitable_agg[cell_index] <- results_annual$months_suitable
annual_relative_abundance_agg[cell_index] <- results_annual$relative_abundance
annual_relative_abundance_agg <- annual_relative_abundance_agg / global(annual_relative_abundance_agg, max, na.rm = TRUE)[[1]]
# values(months_suitable_agg) <- results_annual$months_suitable

library(tidyterra)
ggplot() +
  geom_spatraster(
    data = months_suitable_agg,
    aes(fill = `months suitable`)
  ) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    na.value = grey(0.9)
  ) +
  guides(fill=guide_legend(title = "Months")) +
  theme_minimal() +
  ggtitle("Predicted number of months per year suitable for\nAn. stephensi persistence in microclimate")

ggsave("figures/mechanistic_persistence.png",
       bg = "white",
       width = 9,
       height = 7)

ggplot() +
  geom_spatraster(
    data = annual_relative_abundance_agg,
    aes(fill = `relative abundance (annual)`)
  ) +
  scale_fill_distiller(
    palette = "Greens",
    direction = 1,
    na.value = grey(0.9), 
    trans = "sqrt"
  ) +
  guides(
    fill=guide_legend(title = "Suitability")
  ) +
  theme_minimal() +
  ggtitle("Predicted climatic suitability for An. stephensi\ninside a microclimate")

ggsave("figures/mechanistic_suitability.png",
       bg = "white",
       width = 9,
       height = 7)

# to do:

# run this for An gambiae and for An stephensi

# output rasters with monthly estimates
