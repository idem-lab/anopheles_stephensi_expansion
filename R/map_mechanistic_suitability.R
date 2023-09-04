# Map suitability for An. stephensi and An. gambiae based on Nichemapr
# microclimate model, a cone model of larval habitat, estimated lifehistory
# relationships, and synoptic annual climatic conditions

library(NicheMapR)
library(tidyverse)
library(tidyverse)
library(terra)
library(geodata)
library(tidyterra)

source("R/functions/micro_functions.R")

# load climate-dependent lifehistory functions for An. stephensi and An. gambiae
lifehistory_functions_As <- get_lifehistory_functions("An. stephensi")
lifehistory_functions_Ag <- get_lifehistory_functions("An. gambiae")

# find all the centroids of the climate raster, for the Afro and EMRO regions,
# so we can evaluate this all over the raster

# find and load the global_climate raster downloaded by NicheMapr
gcfolder <- paste(.libPaths()[1], "/gcfolder.rda", sep = "")
load(gcfolder)
nichemapr_climate_raster_mask <- terra::rast(paste0(folder, "/global_climate.nc"))[[1]] * 0 
names(nichemapr_climate_raster_mask) <- "mask"
# crop down to EMRO and AFRO regions

# get a shapefile of the region of interest

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
region_raster_mask_agg <- aggregate(region_raster_mask, 1)
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

# run first for An stephensi microclimate, with permanent water
time_agg <- system.time(
  results_list_As <- sfLapply(
    coords_list,
    calculate_suitability,
    lifehistory_functions = lifehistory_functions_As,
    microclimates = "habitat",
    larval_habitats = "permanent"
    )
)

# time for unaggregated in parallel, ~9h 
time_agg["elapsed"] / (60 * 60)


# then run for An gambiae, non-microclimate, with ephemeral water
time_agg <- system.time(
  results_list_Ag <- sfLapply(
    coords_list,
    calculate_suitability,
    lifehistory_functions = lifehistory_functions_Ag,
    microclimates = "outside",
    larval_habitats = "ephemeral"
  )
)

# time in hours
time_agg["elapsed"] / (60 * 60)

sfStop()


# add on the index to the raster
valid_cell_index <- cellFromXY(region_raster_mask_agg, coords_agg[valid_cells, ])
index_list <- lapply(valid_cell_index, function(x) tibble(cell_index = x))

results_As <- mapply(bind_cols,
                     index_list,
                     results_list_As,
                     SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  `rownames<-`(NULL) %>%
  select(
    cell_index,
    month,
    relative_abundance
  ) %>%
  # scale relative abundance to have maximum value 1
  mutate(
    relative_abundance = relative_abundance / max(relative_abundance, na.rm = TRUE)
  )

results_Ag <- mapply(bind_cols,
                     index_list,
                     results_list_Ag,
                     SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  `rownames<-`(NULL) %>%
  select(
    cell_index,
    month,
    relative_abundance
  ) %>%
  # scale relative abundance to have maximum value 1
  mutate(
    relative_abundance = relative_abundance / max(relative_abundance, na.rm = TRUE)
  )

# work out the index to insert the values
index_raw <- expand.grid(
  cell_index = seq_len(ncell(region_raster_mask_agg)),
  month = 1:12
)

index_As <- index_raw %>%
  left_join(
    results_As,
    by = join_by(cell_index, month)
  ) %>%
  arrange(
    month,
    cell_index,
  )

index_Ag <- index_raw %>%
  left_join(
    results_Ag,
    by = join_by(cell_index, month)
  ) %>%
  arrange(
    month,
    cell_index,
  )


# make an empty monthly raster in which to put results
monthly_relative_abundance_agg <- replicate(12, region_raster_mask_agg, simplify = FALSE) %>%
  do.call(c, .)
names(monthly_relative_abundance_agg) <- month.name

# copy for each species and fill with relevant values
monthly_relative_abundance_agg_As <- monthly_relative_abundance_agg
monthly_relative_abundance_agg_Ag <- monthly_relative_abundance_agg
values(monthly_relative_abundance_agg_As) <- index_As$relative_abundance
values(monthly_relative_abundance_agg_Ag) <- index_Ag$relative_abundance

# fill in NAs for unsuitable areas for An gambiae
monthly_relative_abundance_agg_Ag[is.na(monthly_relative_abundance_agg_Ag)] <- 0
monthly_relative_abundance_agg_Ag <- mask(monthly_relative_abundance_agg_Ag, region_raster_mask_agg)

# save rasters
writeRaster(monthly_relative_abundance_agg_As,
            file = "output/An_stephensi_mechanistic_abundance.tif",
            overwrite = TRUE)

# save rasters
writeRaster(monthly_relative_abundance_agg_Ag,
            file = "output/An_gambiae_mechanistic_abundance.tif",
            overwrite = TRUE)

# plot monthly for An. stephensi
ggplot() +
  geom_spatraster(
    data = monthly_relative_abundance_agg_As,
  ) +
  facet_wrap(~lyr,
             ncol = 3,
             nrow = 4) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = grey(0.9)
  ) +
  labs(fill = "Climatic suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted abundance of An. stephensi per larval habitat")

ggsave("figures/An_stephensi_monthly_mechanistic_suitability.png",
       bg = "white",
       width = 7,
       height = 7)

# plot monthly for An. gambiae
ggplot() +
  geom_spatraster(
    data = monthly_relative_abundance_agg_Ag,
  ) +
  facet_wrap(~lyr,
             ncol = 3,
             nrow = 4) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = grey(0.9)
  ) +
  labs(fill = "Climatic suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted abundance of An. gambiae per larval habitat")

ggsave("figures/An_gambiae_monthly_mechanistic_suitability.png",
       bg = "white",
       width = 7,
       height = 7)

# combine all years and plot
annual_relative_abundance_agg_As <- mean(monthly_relative_abundance_agg_As)
annual_relative_abundance_agg_Ag <- mean(monthly_relative_abundance_agg_Ag)

# plot annual summary for An. stephensi
ggplot() +
  geom_spatraster(
    data = annual_relative_abundance_agg_As,
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = grey(0.9)
  ) +
  labs(fill = "Climatic suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted abundance of An. stephensi per larval habitat")

ggsave("figures/An_stephensi_mechanistic_suitability.png",
       bg = "white",
       width = 7,
       height = 5)

# and for An. gambiae
ggplot() +
  geom_spatraster(
    data = annual_relative_abundance_agg_Ag,
  ) +
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    na.value = grey(0.9)
  ) +
  labs(fill = "Climatic suitability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Predicted abundance of An. gambiae per larval habitat")

ggsave("figures/An_gambiae_mechanistic_suitability.png",
       bg = "white",
       width = 7,
       height = 5)

