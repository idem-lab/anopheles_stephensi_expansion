# Simulating microclimate conditions in An. stephensi habitat over time

# NG's notes for installing NicheMapR (and its dependency gfortran) on Apple silicon
# 
# - install gfortran for macos (apple silicon chips) from here: https://github.com/fxcoudert/gfortran-for-macOS/releases
# - make sure R can see the right gfortran by creating a file ~/.R/Makevars containing the following:
#   
#   FC = /opt/homebrew/bin/gfortran
#   F77 = /opt/homebrew/bin/gfortran
#   FLIBS = -L/opt/homebrew/lib
#   
# -in a fresh R session, do: remotes::install_github("mrke/NicheMapR")

# reinstall NicheMapR from github after 18 August as Mike has fixed a bug when
# running for dry places. To get the specific version I am using:
# remotes::install_github("goldingn/NicheMapR@patch-rainy-NAs")

library(NicheMapR)
library(tidyverse)
library(mgcv)

source("R/functions/micro_functions.R")

# model temperature conditions inside a concrete water container (fully shaded,
# permanent water, 'rock' surface), and compare them to an
# aboveground unshaded situation in the same location

# # height above the ground (metres) at which to calculate conditions, substrate
# # type, and degree of shade in this place - set to 5cm, rock (ie.
# # concrete) and 100% shade to represent a mosquito resting against side of a
# # concrete water tank
# height_m <- 0.05
# soiltype <- 1
# shade_perc <- 100
# wetness_perc <- 100
# 
# # set lat longs of the location
# 
# # placename <- "Awash, Ethiopia"
# # loc <- c(40.142981, 8.9972474)
# 
# # placename <- "Niamey, Niger"
# # loc <- c(2.0840132, 13.5127664)
# 
# # placename <- "Lahore, Pakistan"
# # loc <- c(74.0047326, 31.4831037)
# 
# # placename <- "Bangui, CAR"
# # loc <- c(18.561247, 4.365849)
# 
# # placename <- "Naypyidaw, Myanmar"
# # loc <- c(95.9325102, 19.7470939)
# 
# # placename <- "Iran"
# # loc <- c(57.514125259779256, 28.264071517435397)
# 
placename <- "Salem, India"
loc <- c(78.14653260343869, 11.66683664584138)

# micro <- micro_global(
#   # place
#   loc = loc,
#   timeinterval = 365,
#   # microclimate characteristics
#   Usrhyt = height_m,
#   maxshade = shade_perc,
#   soiltype = soiltype,
#   PCTWET = wetness_perc,
#   runmoist = 0,
#   # make rainfall less mad
#   rainfrac = 0,
#   evenrain = 1
# )
# 
# # plottable dates
# microclimate_temperature_all <- micro$shadmet[, "TALOC"] 
# microclimate_humidity_all <- micro$shadmet[, "RHLOC"] 
# outside_temperature_all <- micro$metout[, "TAREF"]
# outside_humidity_all <- micro$metout[, "RH"]
# 
# # add this in to population model
# microclimate_water_temperature_all <- micro$shadsoil[, "D0cm"]
# 
# temp_range <- range(c(microclimate_temperature_all, outside_temperature_all))
# rh_range <- range(c(microclimate_humidity_all, outside_humidity_all))
# 
# par(mfrow = c(2, 2),
#     oma = c(0, 0, 3, 0))
# 
# # annual profile
# # thin to 2 datapoints per day
# keep <- (micro$dates %% 1) %in% c(0, 0.25, 0.5, 0.75)
# dates <- lubridate::date_decimal(micro$dates[keep] / 365)
# microclimate_temperature <- microclimate_temperature_all[keep]
# microclimate_humidity <- microclimate_humidity_all[keep]
# outside_temperature <- outside_temperature_all[keep]
# outside_humidity <- outside_humidity_all[keep]
# 
# 
# date_range <- range(dates)
# 
# plot(outside_temperature ~ dates,
#      col = "orange",
#      type = "l",
#      lwd = 0.2,
#      ylab = "Temperature (C)",
#      xlab = "",
#      ylim = temp_range,
#      xlim = date_range)
# lines(microclimate_temperature ~ dates, col = "red", lwd = 0.2)
# title("Annual temperature profile")
# 
# 
# plot(outside_humidity ~ dates,
#      col = "lightblue",
#      type = "l",
#      lwd = 0.2,
#      ylab = "Relative humidity (%)",
#      xlab = "",
#      ylim = rh_range,
#      xlim = date_range)
# lines(microclimate_humidity ~ dates,
#       col = "blue",
#       lwd = 0.2)
# title("Annual humidity profile")
# 
# # zoom in on one month
# keep <- seq_along(micro$dates)
# dates <- lubridate::date_decimal(micro$dates[keep] / 365)
# microclimate_temperature <- microclimate_temperature_all[keep]
# microclimate_humidity <- microclimate_humidity_all[keep]
# outside_temperature <- outside_temperature_all[keep]
# outside_humidity <- outside_humidity_all[keep]
# 
# date_range <- lubridate::date_decimal(c(6.5, 6.6) / 12)
# 
# plot(outside_temperature ~ dates,
#      col = "orange",
#      type = "l",
#      lwd = 1,
#      ylab = "Temperature (C)",
#      xlab = "",
#      ylim = temp_range,
#      xlim = date_range)
# lines(microclimate_temperature ~ dates,
#       col = "red",
#       lwd = 2)
# title("Mid-June temperature profile")
# 
# plot(outside_humidity ~ dates,
#      col = "lightblue",
#      type = "l",
#      lwd = 1,
#      ylab = "Relative humidity (%)",
#      xlab = "",
#      ylim = rh_range,
#      xlim = date_range)
# lines(microclimate_humidity ~ dates,
#       col = "blue",
#       lwd = 2)
# title("Mid-June humidity profile")
# 
# title(main = paste("Microclimate conditions in", placename,
#       "
#       either inside a concrete water tank (red, blue) or outside (orange, light blue)"),
#       outer = TRUE)



# port the microclimate conditions into the adult survival model

# load the life history trait functions


# Build a simple two stage model (aquatic stages and adults), with effects of
# density dependence (daily aquatic survival; DAS), water temperature (DAS and
# aquatic development rate; MDR), air temperature (adult survival; DS, and egg
# laying; EFD), and humidity (DS). Construct as a dynamic matrix.

# We have the adult humidity and survival combinations from Bayoh, and we can
# get the other parameters from a combination of studies used in Villena et al.
# and Evans et al. See estimate_lifehistory_functions.R

storage_path <- "data/life_history_params/dehydrated"

# daily adult survival for either An. gambiae or An. stephensi
ds_temp_humid <- rehydrate_lifehistory_function(
  file.path(storage_path, "ds_temp_humid.RDS")
) 

# define An. stephensi only function
ds_function <- function(temperature, humidity) {
  ds_temp_humid(temperature, humidity, species = "An. stephensi")
}

# development rate of aquatic stages as a function of water temperature
mdr_function <- rehydrate_lifehistory_function(
  file.path(storage_path, "mdr_temp_As.RDS")
)

# daily survival probability of aquatic stages as a function of water
# temperature and density of aquatic stages
das_function <- rehydrate_lifehistory_function(
  file.path(storage_path, "das_temp_dens_As.RDS")
)

# daily egg laying as a function of air temperature
efd_function <- rehydrate_lifehistory_function(
  file.path(storage_path, "efd_temp_As.RDS")
)


# # get hourly adult survival probabilities
# period_days <- 1 / 24
# survival_microclimate_all <- ds_function(temperature = microclimate_temperature_all,
#                                          humidity = microclimate_humidity_all) ^ period_days
# 
# survival_outside_all <- ds_function(temperature = outside_temperature_all,
#                                     humidity = outside_humidity_all) ^ period_days
# 
# survival_microclimate <- survival_microclimate_all[keep]
# survival_outside <- survival_outside_all[keep]
# par(mfrow = c(1, 1))
# plot(survival_outside ~ dates,
#      col = "light green",
#      type = "l",
#      lwd = 2,
#      ylab = "Hourly probability of survival",
#      xlab = "",
#      ylim = range(c(survival_microclimate, survival_outside)),
#      xlim = date_range)
# lines(survival_microclimate ~ dates,
#       col = "dark green",
#       lwd = 2)



# das_function is a function of temperature and density (larvae per 250ml
# water), the others are functions of temperature

# given our microclimate data, we can now compute these parameters

# For the larval stages (mdr, das) we can use water temperature. The
# experiements use air temperature, but small volumes of water and high humidity
# so that they track with the air temperatures, but in our microclimate there
# will be a buffering effect

conditions <- model_climatic_conditions(loc)
scaling <- 1/24
mdr <- mdr_function(conditions$habitat$water_temperature) ^ scaling
das <- das_function(conditions$habitat$water_temperature, density = 0) ^ scaling
das64 <- das_function(conditions$habitat$water_temperature, density = 64) ^ scaling
efd <- efd_function(conditions$habitat$air_temperature) ^ scaling
ds <- ds_function(temperature = conditions$habitat$air_temperature,
                  humidity = conditions$habitat$humidity) ^ scaling


par(mfrow = c(2, 2),
    oma = rep(0, 4),
    mar = c(3, 3, 4, 1) + 0.1)
plot(ds ~ conditions$habitat$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "", ylab = "",
     main = "adult survival prob")
plot(das ~ conditions$habitat$day,
     type = "l",
     xlim = c(10, 14),
     ylim = range(das, das64),
     xlab = "", ylab = "",
     main = "aquatic survival prob \n(low and high density)")
lines(das64 ~ conditions$habitat$day,
      lty = 2)
plot(mdr ~ conditions$habitat$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "day of year", ylab = "",
     main = "larval development rate")
plot(efd ~ conditions$habitat$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "day of year", ylab = "",
     main = "egg laying rate")

# put this into a population dynamic simulation model

conditions <- model_climatic_conditions(loc)
states <- simulate_population(conditions$habitat)

par(mfrow = c(2, 1),
    mar = c(4, 4, 1, 2) + 0.1)
plot(states[, 1] ~ conditions$habitat$day,
     ylab = "aquatic stages",
     xlab = "",
     ylim = c(0, max(states[, 1])),
     type = "l")
plot(states[, 2] ~ conditions$habitat$day,
     ylab = "adults",
     type = "l",
     ylim = c(0, max(states[, 2])),
     xlab = "day")

# placename <- "Awash, Ethiopia"
# loc <- c(40.142981, 8.9972474)
suit_awash <- calculate_stephensi_suitability(c(40.142981, 8.9972474))

# placename <- "Niamey, Niger"
# loc <- c(2.0840132, 13.5127664)
suit_niamey <- calculate_stephensi_suitability(c(2.0840132, 13.5127664))

# placename <- "Bangui, CAR"
# loc <- c(18.561247, 4.365849)
suit_bangui <- calculate_stephensi_suitability(c(18.561247, 4.365849))

suit <- bind_rows(
  mutate(suit_awash, location = "Awash, Ethiopia"),
  mutate(suit_niamey, location = "Niamey, Niger"),
  mutate(suit_bangui, location = "Bangui, CAR"),
)

library(ggplot2)

suit %>%
  mutate(
    month = factor(month, levels = unique(month))
  ) %>%
  ggplot(
    aes(x = month,
        y = relative_abundance,
        group = microclimate,
        colour = microclimate)
  ) +
  geom_line() +
  facet_grid(larval_habitat ~ location) +
  theme_minimal()


suit %>%
  mutate(
    month = factor(month, levels = unique(month))
  ) %>%
  ggplot(
    aes(x = month,
        y = persistence,
        group = microclimate,
        colour = microclimate)
  ) +
  geom_line() +
  facet_grid(larval_habitat ~ location) +
  theme_minimal()


suit %>%
  filter(
    microclimate == "habitat"
  ) %>%
  mutate(
    month = factor(month, levels = unique(month))
  ) %>%
  ggplot(
    aes(
      x = month,
      y = relative_abundance,
      group = larval_habitat,
      colour = larval_habitat
    )
  ) +
  geom_line() +
  facet_grid(~ location) +
  theme_minimal()


# run this over all pixels to map year-round microclimate (and outside
# microclimate) suitability





# find all the centroids of the climate raster, for the Afro and EMRO regions,
# so we can evaluate this all over the raster

# find and load the global_climate raster
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


# aggregate it 2x to get a lower resolution set of coordinates for second pass -
# original resolution is ~18.5km at the equator
region_raster_mask_agg <- aggregate(region_raster_mask, 4)
coords_agg <- xyFromCell(region_raster_mask_agg, cells(region_raster_mask_agg))
nrow(coords_agg)


# find and skip NAs in the global climate raster at these locations (use raster,
# as terra gives slightly different locations?)
nichemapr_climate_raster <- raster::brick(paste0(folder, "/global_climate.nc"))
vals <- raster::extract(nichemapr_climate_raster, coords_agg)
fine <- apply(!is.na(vals), 1, all)
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
sfExport(list = list("model_climatic_conditions",
                     "simulate_population",
                     "mdr_function",
                     "das_function",
                     "efd_function",
                     "ds_function",
                     "ds_temp_humid",
                     "iterate_state",
                     "create_matrix",
                     "summarise_dynamics"))
time_agg <- system.time(
  results_list <- sfLapply(coords_list, calculate_stephensi_suitability)
)



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
    microclimate == "habitat"
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

# note: we could consider switching to the terraclimate inputs and run over
# actual time? It might not be much slower to run for a decade, given overheads
# in the microclimate model

# need to check model with Mike Kearney, and consider doing an outside water
# source model too

# Increase the water volume size when determining the persistence suitability,
# because populations will consist of multiple water tanks

# # try recreating the timeseries of monthly abundance in Whittaker et al. (2023)
# # https://doi.org/10.1073/pnas.2216142120
# jalalabad_loc <- rev(c(34.4198911, 70.4303445))
# lahore_loc <- rev(c(31.4831037, 74.0047326))
# djibouti_loc <- rev(c(11.5747928, 43.0778374))
# bandar_abbas_loc <- rev(c(27.1973499, 55.9655767))
# naypyidaw_loc <- rev(c(19.7470939, 95.9325102))
# 
# afghanistan <- calculate_stephensi_suitability(jalalabad_loc) %>%
#   mutate(location = "Afghanistan")
# pakistan <- calculate_stephensi_suitability(lahore_loc) %>%
#   mutate(location = "Pakistan")
# djibouti <- calculate_stephensi_suitability(djibouti_loc) %>%
#   mutate(location = "Djibouti")
# iran <- calculate_stephensi_suitability(bandar_abbas_loc) %>%
#   mutate(location = "Iran")
# myanmar <- calculate_stephensi_suitability(naypyidaw_loc) %>%
#   mutate(location = "Myanmar")
# 
# month_letter <- substr(month.abb, 1, 1)
# 
# bind_rows(
#   afghanistan,
#   pakistan,
#   djibouti,
#   iran,
#   myanmar
# ) %>%
#   filter(
#     microclimate == "habitat"
#   ) %>%
#   mutate(
#     relative_abundance = relative_abundance / max(relative_abundance),
#     month = factor(month.abb[month],
#                    levels = month.abb)
#   ) %>%
#   ggplot(
#     aes(
#       x = month,
#       y = relative_abundance,
#       group = location
#     )
#   ) +
#   scale_x_discrete(
#     labels = month_letter
#   ) +
#   geom_line() +
#   facet_grid(larval_habitat ~ location) +
#   theme_minimal()

# read in Whittaker et al paper information
whittaker_papers <- read_csv("https://raw.githubusercontent.com/goldingn/stephenseasonality/main/data/systematic_review_results/extracted_entomological_data.csv") %>%
  select(`New ID`,
         id = `Time Series ID`,
         Year,
         Author,
         Title,
         Country,
         `Admin 1`,
         `Admin 2`) %>%
  distinct()

# read in Whittaker et al extracted data (from forked repo to safeguard against changes)
whittaker_data <- get_rds("https://github.com/goldingn/stephenseasonality/raw/main/data/systematic_review_results/metadata_and_processed_unsmoothed_counts.rds")
whittaker_admin1 <- get_rds("https://github.com/goldingn/stephenseasonality/raw/main/data/admin_units/simplified_admin1.rds")
whittaker_admin2 <- get_rds("https://github.com/goldingn/stephenseasonality/raw/main/data/admin_units/simplified_admin2.rds")

# punjab appears twice here
admin1 <- whittaker_admin1 %>%
  select(
    country = NAME_0,
    admin1 = NAME_1,
    shape_admin1 = geometry
  )

admin2 <- whittaker_admin2 %>%
  select(
    country = NAME_0,
    admin1 = NAME_1,
    admin2 = NAME_2,
    shape_admin2 = geometry
  )
library(sf)
  

# get centroids for regions

whittaker_tidied <- whittaker_data %>%
  as_tibble() %>%
  # only keep complete annual timeseries
  rowwise() %>%
  filter(
    !any(across(
      any_of(month.abb),
      ~is.na(.x)
    ))
  ) %>%
  # only keep timeseries with at least 100 mossies per year (to remove noise)
  rowwise() %>%
  mutate(
    n_years = end - start + 1,
    total = sum(across(any_of(month.abb)), na.rm = TRUE),
    total_per_year = total / n_years,
    .after = id
  ) %>%
  filter(
    # want them to average at least 25 per month (it would ideally be more, but
    # trying to keep some data in :/ )
    total >= 300
  ) %>%
  # keep only places with multiple years of data
  rowwise() %>%
  mutate(
    years = list(seq(start, end))
  ) %>%
  group_by(country, admin1, admin2, city) %>%
  mutate(
    location_n_years = n_distinct(unlist(years)),
    .after = id
  ) %>%
  ungroup() %>%
  filter(
    location_n_years >= 3
  ) %>%
  # keep only timeseries in places with at least 2 timeseries (to understand
  # variability)
  group_by(country, admin1, admin2, city) %>%
  mutate(
    timeseries = n()
  ) %>%
  ungroup() %>%
  filter(
    timeseries >= 2
  ) %>%
  # join on spatial data and find centroids
  left_join(
    admin1,
    by = c("country", "admin1")
  ) %>%
  left_join(
    admin2,
    by = c("country", "admin1", "admin2")
  ) %>%
  mutate(
    shape_admin = case_when(
      is.na(admin2) ~ shape_admin1,
      .default = shape_admin2
    )
  ) %>%
  mutate(
    coords = st_centroid(shape_admin)
  ) %>%
  mutate(
    placename = case_when(
      is.na(admin2) ~ paste(admin1, country, sep = ", "),
      .default = paste(admin2, country, sep = ", ")
    ),
    .after = id
  ) %>%
  select(
    -shape_admin1,
    -shape_admin2,
    -shape_admin
  )

# unique locations with the most data
best_data <- whittaker_tidied %>%
  group_by(placename, city, coords) %>%
  summarise(
    timeseries = timeseries[1],
    mossies_per_timeseries = sum(total) / timeseries,
    .groups = "drop"
  ) %>%
  arrange(desc(mossies_per_timeseries))

whittaker_subset <- whittaker_tidied %>%
  filter(placename %in% best_data$placename)

whittaker_coords_list <- best_data %>%
  pull(coords) %>%
  st_coordinates() %>%
  as_tibble() %>%
  split(seq_len(nrow(.)))

whittaker_results_list <- lapply(whittaker_coords_list, calculate_stephensi_suitability)

index_list <- lapply(best_data$placename, function(x) tibble(placename = x))

whittaker_results <- mapply(bind_cols, index_list, whittaker_results_list, SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  `rownames<-`(NULL)

whittaker_obs_pred <- whittaker_subset %>%
  pivot_longer(
    cols = month.abb,
    names_to = "month",
    values_to = "abundance"
  ) %>%
  mutate(
    month_id = match(month, month.abb)
  ) %>%
  left_join(
    filter(
      whittaker_results,
      microclimate == "habitat",
      larval_habitat == "permanent"
    ) %>%
      select(-microclimate,
             -larval_habitat),
    by = c("placename", month_id = "month")
  ) %>%
  rename(
    relative_abundance_habitat_permanent = relative_abundance
  ) %>%
  left_join(
    filter(
      whittaker_results,
      microclimate == "habitat",
      larval_habitat == "ephemeral"
    ) %>%
      select(-microclimate,
             -larval_habitat),
    by = c("placename", month_id = "month")
  ) %>%
  rename(
    relative_abundance_habitat_ephemeral = relative_abundance
  ) %>%
  group_by(id) %>%
  mutate(
    abundance = abundance / mean(abundance, na.rm = TRUE),
    across(
      starts_with("relative_abundance"),
      ~ .x / mean(.x)
    )
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = starts_with("relative_abundance"),
    names_to = "Modelled",
    names_prefix = "relative_abundance_habitat_",
    values_to = "relative_abundance"
  ) %>%
  mutate(
    Modelled = case_when(
      Modelled == "permanent" ~ "Permanent water",
      Modelled == "ephemeral" ~ "Ephemeral water"
    ),
    # make permanent a solid line, and ephemeral a dashed line
    Modelled = factor(
      Modelled,
      levels = c("Permanent water",
                 "Ephemeral water")
    )
  ) %>%
  # tidy up the placename for Delhi
  mutate(
    placename = gsub("NCT of Delhi",
                     "Delhi",
                     placename)
  ) %>%
  # add the year range for the location
  group_by(placename) %>%
  mutate(
    year_range = paste(min(start), max(end), sep = "-")
  ) %>%
  ungroup() %>%
  mutate(
    panel_name = sprintf("%s %s\n(%s)",
                         city,
                         placename, 
                         year_range),
    panel_name = factor(panel_name,
                        levels = rev(unique(panel_name)))
  ) %>%
  # rename the urbanness
  rename(
    `Data` = city
  )

# get the peak rain periods to plot
rain_months <- lapply(whittaker_coords_list, get_rain_months)
index_list <- lapply(best_data$placename, function(x) tibble(placename = x))

rain_months_plot <- mapply(bind_cols, index_list, rain_months, SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  `rownames<-`(NULL) %>%
  group_by(placename) %>%
  mutate(
    peak = rainfall == max(rainfall),
    `Proportion rainfall` = rainfall / sum(rainfall),
    near_peak = `Proportion rainfall` >= 0.1,
    abundance = NA,
    id = NA,
    Data = NA,
    placename = gsub("NCT of Delhi", "Delhi", placename)
  ) %>%
  left_join(
    whittaker_obs_pred %>%
      select(placename, panel_name) %>%
      distinct(),
    by = "placename"
  )

plot <- whittaker_obs_pred %>%
  mutate(
    month = factor(month, levels = month.abb)
  ) %>%
  ggplot(
    aes(
      y = abundance,
      x = month,
      group = id,
      colour = Data
    )
  ) +
  geom_vline(
    data = rain_months_plot,
    aes(
      xintercept = month,
      alpha = `Proportion rainfall`
    ),
    colour = "skyblue",
    linewidth = 7
  ) +
  scale_alpha_continuous(labels = scales::percent,
                         trans = "exp") +
  geom_line(alpha = 0.3,
            linewidth = 0.5) +
  geom_point(
    alpha = 0.3,
    # pch = 15,
    size = 2
  ) +
  scale_color_brewer(palette = "Set1") +
  geom_line(
    aes(
      y = relative_abundance,
      group = Modelled,
      linetype = Modelled
    ),
    colour = "black",
    linewidth = 1
  ) +
  facet_wrap(~panel_name) +
  ylab("Relative abundance") +
  xlab("") +
  theme_minimal() +
  theme(
    axis.line = element_line(color='black'),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

ggsave(
  filename = "figures/whittaker_comparison.png",
  plot = plot,
  width = 8,
  height = 6,
  bg = "white"
)

# format the papers to include in the slide
whittaker_papers %>%
  filter(
    id %in% unique(whittaker_obs_pred$id)
  ) %>%
  mutate(
    Place = case_when(
      is.na(`Admin 2`) ~ `Admin 1`,
      .default = `Admin 2`
    ),
    Place = gsub("NCT of Delhi", "Delhi", Place),
    Place = paste(Place, Country, sep = ", "),
    Paper = sprintf("%s (%s)", Author, Year)
  ) %>%
  select(Place, Paper) %>%
  distinct() %>%
  group_by(Place) %>%
  summarise(
    text = sprintf("%s:\t%s",
                   Place,
                   paste(unlist(Paper), collapse = "; ")
    ),
    .groups = "drop"
  ) %>%
  distinct() %>%
  pull(text) %>%
  paste(collapse = "\n") %>%
  write_file(file = "figures/whittaker_comparison_text.txt")
   


awash_loc <- c(40.142981, 8.9972474)
addis_loc <- c(38.75960099237313, 8.982497645041477)
# loc <- addis_loc
loc <- addis_loc
climate <- format_climatic_data(loc)

variables_keep <- c("Air temperature (C)",
                    "Water temperature (C)",
                    "Humidity (%)",
                    "Rainfall (mm)",
                    "Pooled water\n(relative area)")
variables_col <- RColorBrewer::brewer.pal(9, "Set1")[c(5, 8, 3, 2, 2)]


# plot the microclimate hourly and annually
plot_micro_hourly <- climate %>%
  # plot only these variables, and in the right order
  filter(
    variable %in% variables_keep[c(1:3)]
  ) %>%
  # one day either side of presentation day, plus 3h time difference from GMT
  filter(
    date >= (as.Date("2023-09-19") - 1),
    date <= (as.Date("2023-09-19") + 1),
  ) %>%
  mutate(
    variable = factor(variable,
                      levels = variables_keep)
  ) %>%
  # express rainfall as mm/day
  mutate(
    multiplier = case_when(
      variable == "Rainfall (mm)" ~ 24,
      .default = 1
    ),
    across(
      starts_with("value"),
      ~ .x * multiplier
    ),
    # convert to hours since start
    hour = round((day - min(day)) * 24)
  ) %>%
  ggplot(
    aes(
      x = hour,
      y = value,
      linetype = which,
      colour = variable,
      fill = variable
    )
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_colour_manual(
    labels = variables_keep,
    values = variables_col,
    guide = "none"
  ) +
  scale_fill_manual(
    labels = variables_keep,
    values = variables_col,
    guide = "none"
  ) +
  facet_grid(rows = "variable", 
    # ~ variable ~, ~  which,
    scales = "free",
    switch = "y"
  ) +
  ylab("") +
  scale_x_continuous(
    breaks = seq(0, 3 * 24, by = 6),
    # label has hour of the day, and shift to Ethiopia time
    labels = function(hours) {
      (hours - 1) %% 24 + 1
    }
  ) +
  scale_linetype(name = "") +
  theme_minimal() +
  theme(
    strip.placement = "outside"
  )

plot_micro_year <- climate %>%
  # plot only these variables, and in the right order
  filter(
    variable %in% variables_keep,
    which == "microclimate"
  ) %>%
  mutate(
    variable = factor(variable,
                      levels = variables_keep)
  ) %>%
  # summarise by week for plotting
  group_by(which, week, variable) %>%
  summarise(
    date = mean(date),
    value_mean = mean(value),
    value_upper = quantile(value, 0.9),
    value_lower = quantile(value, 0.1),
    .groups = "drop"
  ) %>%
  # express rainfall as mm/day
  mutate(
    multiplier = case_when(
      variable == "Rainfall (mm)" ~ 24,
      .default = 1
    ),
    across(
      starts_with("value"),
      ~ .x * multiplier
    )
  ) %>%
  ggplot(
    aes(
      x = date,
      y = value_mean,
      ymax = value_upper,
      ymin = value_lower,
      # linetype = which,
      colour = variable,
      fill = variable
    )
  ) +
  geom_ribbon(
    alpha = 0.2,
    linewidth = 0.1
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_colour_manual(
    labels = variables_keep,
    values = variables_col,
    guide = "none"
  ) +
  scale_fill_manual(
    labels = variables_keep,
    values = variables_col,
    guide = "none"
  ) +
  scale_x_date(
    date_labels = "%b"
  ) +
  facet_grid(
    rows = "variable",
    scales = "free",
    switch = "y"
  ) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(
    strip.placement = "outside"
  )

ggsave("figures/microclimate_hourly.png",
       plot_micro_hourly,
       bg = "white",
       height = 5,
       width = 5)

ggsave("figures/microclimate_year.png",
       plot_micro_year,
       bg = "white",
       height = 5,
       width = 5)


# tidy up visualisation of climate and lifehistory timeseries (ggplot code from
# modelling conditions and larval habitat on different timeframes)

# clean up script

# remake rasters including larval habitat

