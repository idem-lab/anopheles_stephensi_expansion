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

# # Note it's possible to use ERA5 grided past weather data, but it's slooow
# # for era5: install R package needed for linking to ERA5 data
# remotes::install_github("dklinges9/mcera5")
# # folder and file prefix with local ERA5 data extracted via the mcera5 package (no trailing forward slash)
# era5_dir <- "~/build/era5_data"
# era5_prefix <- "era5"
# era5_dir_prefix <- file.path(era5_dir, era5_prefix)
# if (!dir.exists(era5_dir)) {
#   dir.create(era5_dir)
# }
# # set copernicus API key
# uid <- "188032"
# cds_api_key <- "these can be found at the bottom of the user profile"
# # after registering here: https://cds.climate.copernicus.eu/user/register
# ecmwfr::wf_set_key(user = uid, key = cds_api_key, service = "cds")
# # you also need to accept the license terms here:
# # https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products
# # define a bounding box around the location for download
# xmin <- 140
# xmax <- 144
# ymin <- -36
# ymax <- -32
# # and time period
# st_time <- "2016-01-01"
# en_time <- "2018-12-31"
# # build a request
# req <- build_era5_request(xmin = xmin,
#                           xmax = xmax,
#                           ymin = ymin,
#                           ymax = ymax,
#                           start_time = st_time,
#                           end_time = en_time,
#                           outfile_name = era5_prefix)
# # download the data
# request_era5(request = req, uid = uid, out_path = era5_dir)

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
      group = location
    )
  ) +
  geom_line() +
  facet_grid(larval_habitat ~ location) +
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

# try recreating the timeseries of monthly abundance in Whittaker et al. (2023)
# https://doi.org/10.1073/pnas.2216142120
jalalabad_loc <- rev(c(34.4198911, 70.4303445))
lahore_loc <- rev(c(31.4831037, 74.0047326))
djibouti_loc <- rev(c(11.5747928, 43.0778374))
bandar_abbas_loc <- rev(c(27.1973499, 55.9655767))
naypyidaw_loc <- rev(c(19.7470939, 95.9325102))

afghanistan <- calculate_stephensi_suitability(jalalabad_loc) %>%
  mutate(location = "Afghanistan")
pakistan <- calculate_stephensi_suitability(lahore_loc) %>%
  mutate(location = "Pakistan")
djibouti <- calculate_stephensi_suitability(djibouti_loc) %>%
  mutate(location = "Djibouti")
iran <- calculate_stephensi_suitability(bandar_abbas_loc) %>%
  mutate(location = "Iran")
myanmar <- calculate_stephensi_suitability(naypyidaw_loc) %>%
  mutate(location = "Myanmar")

month_letter <- substr(month.abb, 1, 1)

bind_rows(
  afghanistan,
  pakistan,
  djibouti,
  iran,
  myanmar
) %>%
  filter(
    microclimate == "habitat"
  ) %>%
  mutate(
    relative_abundance = relative_abundance / max(relative_abundance),
    month = factor(month.abb[month],
                   levels = month.abb)
  ) %>%
  ggplot(
    aes(
      x = month,
      y = relative_abundance,
      group = location
    )
  ) +
  scale_x_discrete(
    labels = month_letter
  ) +
  geom_line() +
  facet_grid(larval_habitat ~ location) +
  theme_minimal()

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

whittaker_papers %>%
  filter(`Admin 2` == "Kheda")

whittaker_data %>% filter(admin2 == "Kheda")

# read in Whittaker et al extracted data (from forked repo to safeguard against changes)
whittaker_data <- get_rds("https://github.com/goldingn/stephenseasonality/raw/main/data/systematic_review_results/metadata_and_processed_unsmoothed_counts.rds")
whittaker_admin1 <- get_rds("https://github.com/goldingn/stephenseasonality/raw/main/data/admin_units/simplified_admin1.rds")
whittaker_admin2 <- get_rds("https://github.com/goldingn/stephenseasonality/raw/main/data/admin_units/simplified_admin2.rds")

# punjab apdasrs twice here
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
    coords = st_centroid(shape_admin),
    placename = paste(admin2, admin1, country, sep = ", ")
  ) %>%
  select(
    -shape_admin1,
    -shape_admin2,
    -shape_admin
  )

# find the locations with the most years (at least 4), and the most mosquitoes
# caught per year
best_data <- whittaker_tidied %>%
  group_by(country, admin1, admin2, city, coords) %>%
  summarise(
    total = sum(across(any_of(month.abb)), na.rm = TRUE),
    years = n(),
    year_average = total / years,
    .groups = "drop") %>%
  filter(
    !is.na(admin2),
    years >= 4
  ) %>%
  arrange(desc(year_average)) %>%
  mutate(
    placename = paste(admin2, admin1, country, sep = ", "),
    .before = everything()
  )

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
  left_join(
    filter(
      whittaker_results,
      microclimate == "outside",
      larval_habitat == "permanent"
    ) %>%
      select(-microclimate,
             -larval_habitat),
    by = c("placename", month_id = "month")
  ) %>%
  rename(
    relative_abundance_outside_permanent = relative_abundance
  ) %>%
  left_join(
    filter(
      whittaker_results,
      microclimate == "outside",
      larval_habitat == "ephemeral"
    ) %>%
      select(-microclimate,
             -larval_habitat),
    by = c("placename", month_id = "month")
  ) %>%
  rename(
    relative_abundance_outside_ephemeral = relative_abundance
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
  # keep only the higher resolution
  filter(!is.na(admin2))

ylims <- whittaker_obs_pred %>%
  select(
    abundance,
    starts_with("relative_abundance")
  ) %>%
  as.matrix() %>%
  range(na.rm = TRUE)

plot <- whittaker_obs_pred %>%
  mutate(
    month = factor(month, levels = month.abb)
  ) %>%
  ggplot(
    aes(
      y = abundance,
      x = month,
      group = id,
      colour = city
    )
  ) +
  geom_line(
  ) +
  geom_point(
  ) +
  geom_line(
    aes(
      y = relative_abundance_habitat_permanent
    ),
    colour = "black",
    linetype = 1
  ) +
  geom_line(
    aes(
      y = relative_abundance_habitat_ephemeral
    ),
    colour = "black",
    linetype = 2
  ) +
  # geom_line(
  #   aes(
  #     y = relative_abundance_outside_permanent
  #   ),
  #   colour = "grey40",
  #   linetype = 1
  # ) +
  # geom_line(
  #   aes(
  #     y = relative_abundance_outside_ephemeral
  #   ),
  #   colour = "grey40",
  #   linetype = 2
  # ) +
  coord_cartesian(
    ylim = ylims
  ) +
  facet_wrap(~placename) +
  theme_minimal()

# run two versions of the population model - one with permanent larval habitat, and one from
# rain-fed larval habitat

# 1. smoothly interpolate rainfall (from monthly data)
# 2a. investigate using rainfall as linear in abundance 

# 2b. implement Owen's puddle model for amount of larval habitat (either
# microclimate or outdoor evaporation metrics)
# 3. preprocess larval habitat area for both models
# 4. run mosquito populations with both sorts of larval habitat (permanent or
# ephemeral)


placename <- "Salem, India"
loc <- c(78.14653260343869, 11.66683664584138)

# solve the cone model forward through time on an hourly timestep
cond <- model_climatic_conditions(loc)
conditions <- cond$outside
# conditions$rainfall <- cond$ephemeral_larval_habitat
# conditions$windspeed <- micro$shadmet[, "VLOC"]
# conditions$altitude <- 0

# The inflow multiplier is the ratio of rainfall catchment to the maximum larval
# surface area. When it is large, the larval habitat is more likely to max-out.
# This can potentially be tweaked to replicate the type of waater body (a
# puddle, where there's likely to be no run-on that can't fill up) versus a
# rainwater tank, where the additional catchment is likely to be significant

larval_habitat_in <- simulate_ephemeral_habitat(conditions = cond$habitat,
                                                inflow_multiplier = 6)
larval_habitat_out <- simulate_ephemeral_habitat(conditions = cond$outside,
                                                 inflow_multiplier = 6)

cone_volume_to_surface(cone_depth_to_volume(0.1))

par(mfrow = c(2, 1))
plot(rainfall ~ day,
     type = "l",
     col = "blue",
     data = cond$habitat)
plot(larval_habitat_in ~ day,
     ylim = range(c(larval_habitat_in, larval_habitat_out)),
     type = "l",
     col = "pink",
     data = cond$habitat)
lines(larval_habitat_out ~ day,
     col = "purple",
     data = cond$habitat)
abline(h = c(larval_habitat_in[1],
             larval_habitat_out[1]),
       lty = 2)

suit <- calculate_stephensi_suitability(loc)

suit %>%
  mutate(
    month = factor(month, levels = unique(month))
  ) %>%
  ggplot(
    aes(x = month,
        y = relative_abundance,
        group = larval_habitat,
        colour = larval_habitat)
  ) +
  geom_line() +
  facet_grid(~ microclimate) +
  theme_minimal()

# to do:

# resolve Niamey error

# tidy up plotting against Whittaker data

# tweak model to fit Whittaker data

# tidy up visualisation of climate and lifehistory timeseries (ggplot code from
# modelling conditions and larval habitat on different timeframes)

# clean up script

# remake rasters including larval habitat
