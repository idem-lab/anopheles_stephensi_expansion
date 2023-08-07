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

# model temperature conditions inside a concrete water container (fully shaded,
# permanent water, 'rock' surface), and compare them to an
# aboveground unshaded situation in the same location

# height above the ground (metres) at which to calculate conditions, substrate
# type, and degree of shade in this place - set to 5cm, rock (ie.
# concrete) and 100% shade to represent a mosquito resting against side of a
# concrete water tank
height_m <- 0.05
soiltype <- 1
shade_perc <- 100
wetness_perc <- 100

# set lat longs of the location

# placename <- "Awash, Ethiopia"
# loc <- c(40.142981, 8.9972474)

placename <- "Niamey, Niger"
loc <- c(2.0840132, 13.5127664)

# placename <- "Bangui, CAR"
# loc <- c(18.561247, 4.365849)


micro <- micro_global(
  # place
  loc = loc,
  timeinterval = 365,
  # microclimate characteristics
  Usrhyt = height_m,
  maxshade = shade_perc,
  soiltype = soiltype,
  PCTWET = wetness_perc,
  runmoist = 0
  # puddle characteristics
  # rainmult = rainmult,
  # maxpool = maxpool,
  # evenrain = 1
)

# plottable dates
microclimate_temperature_all <- micro$shadmet[, "TALOC"] 
microclimate_humidity_all <- micro$shadmet[, "RHLOC"] 
outside_temperature_all <- micro$metout[, "TAREF"]
outside_humidity_all <- micro$metout[, "RH"]

# add this in to population model
microclimate_water_temperature_all <- micro$shadsoil[, "D0cm"]

temp_range <- range(c(microclimate_temperature_all, outside_temperature_all))
rh_range <- range(c(microclimate_humidity_all, outside_humidity_all))

par(mfrow = c(2, 2),
    oma = c(0, 0, 3, 0))

# annual profile
# thin to 2 datapoints per day
keep <- (micro$dates %% 1) %in% c(0, 0.25, 0.5, 0.75)
dates <- lubridate::date_decimal(micro$dates[keep] / 365)
microclimate_temperature <- microclimate_temperature_all[keep]
microclimate_humidity <- microclimate_humidity_all[keep]
outside_temperature <- outside_temperature_all[keep]
outside_humidity <- outside_humidity_all[keep]


date_range <- range(dates)

plot(outside_temperature ~ dates,
     col = "orange",
     type = "l",
     lwd = 0.2,
     ylab = "Temperature (C)",
     xlab = "",
     ylim = temp_range,
     xlim = date_range)
lines(microclimate_temperature ~ dates, col = "red", lwd = 0.2)
title("Annual temperature profile")


plot(outside_humidity ~ dates,
     col = "lightblue",
     type = "l",
     lwd = 0.2,
     ylab = "Relative humidity (%)",
     xlab = "",
     ylim = rh_range,
     xlim = date_range)
lines(microclimate_humidity ~ dates,
      col = "blue",
      lwd = 0.2)
title("Annual humidity profile")

# zoom in on one month
keep <- seq_along(micro$dates)
dates <- lubridate::date_decimal(micro$dates[keep] / 365)
microclimate_temperature <- microclimate_temperature_all[keep]
microclimate_humidity <- microclimate_humidity_all[keep]
outside_temperature <- outside_temperature_all[keep]
outside_humidity <- outside_humidity_all[keep]

date_range <- lubridate::date_decimal(c(6.5, 6.6) / 12)

plot(outside_temperature ~ dates,
     col = "orange",
     type = "l",
     lwd = 1,
     ylab = "Temperature (C)",
     xlab = "",
     ylim = temp_range,
     xlim = date_range)
lines(microclimate_temperature ~ dates,
      col = "red",
      lwd = 2)
title("Mid-June temperature profile")

plot(outside_humidity ~ dates,
     col = "lightblue",
     type = "l",
     lwd = 1,
     ylab = "Relative humidity (%)",
     xlab = "",
     ylim = rh_range,
     xlim = date_range)
lines(microclimate_humidity ~ dates,
      col = "blue",
      lwd = 2)
title("Mid-June humidity profile")

title(main = paste("Microclimate conditions in", placename,
      "
      either inside a concrete water tank (red, blue) or outside (orange, light blue)"),
      outer = TRUE)



# port the microclimate conditions into the adult survival model

# load the life history trait functions


# Build a simple two stage model (aquatic stages and adults), with effects of
# density dependence (daily aquatic survival; DAS), water temperature (DAS and
# aquatic development rate; MDR), air temperature (adult survival; DS, and egg
# laying; EFD), and humidity (DS). Construct as a dynamic matrix.

# We have the adult humidity and survival combinations from Bayoh, and we can
# get the other parameters from a combination of studies used in Villena et al.
# and Evans et al. See estimate_lifehistory_functions.R

# reload these as functions from saved objects
rehydrate_lifehistory_function <- function(path_to_object) {
  object <- readRDS(path_to_object)
  do.call(`function`,
          list(object$arguments,
               body(object$dummy_function)))
}

storage_path <- "data/life_history_params/dehydrated"

# daily adult survival for either An. gambiae or An. stephensi
ds_temp_humid <- rehydrate_lifehistory_function(
  file.path(storage_path, "ds_temp_humid.RDS")
) 

# An. stephensi only function
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

# 
# 
# 
# survival_model <- readRDS("data/survival_model/survival_model.RDS")
# 
# # function to compute adult survival at temperature and humidity combinations
# survival <- function (temperature = 20,
#                       humidity = 100,
#                       use_wild_correction = TRUE,
#                       period_days = 1) {
#   
#   # construct dataframe for prediction
#   df_pred <- data.frame(temperature = temperature,
#                         humidity = humidity,
#                         sex = "F",
#                         off = 0)
#   
#   # predict cumulative survival probability after one day
#   survival <- 1 - predict(survival_model, df_pred, type = "response")
#   
#   if (use_wild_correction) {
#     wild_correction <- attr(survival_model, "wild_correction")
#     survival <- survival * wild_correction
#   }
#   
#   # correct to given time period over which to calculate survival
#   survival <- survival ^ period_days
#   
#   survival
#   
# }

# get hourly adult survival probabilities
period_days <- 1 / 24
survival_microclimate_all <- ds_function(temperature = microclimate_temperature_all,
                                         humidity = microclimate_humidity_all) ^ period_days

survival_outside_all <- ds_function(temperature = outside_temperature_all,
                                    humidity = outside_humidity_all) ^ period_days

survival_microclimate <- survival_microclimate_all[keep]
survival_outside <- survival_outside_all[keep]
par(mfrow = c(1, 1))
plot(survival_outside ~ dates,
     col = "light green",
     type = "l",
     lwd = 2,
     ylab = "Hourly probability of survival",
     xlab = "",
     ylim = range(c(survival_microclimate, survival_outside)),
     xlim = date_range)
lines(survival_microclimate ~ dates,
      col = "dark green",
      lwd = 2)

# given a latitude and a longitude, model the temperuture and humidity profile
# at an hourly resolution over an average year, both inside a hypothetical
# concrete water tank, and out in the open air
model_conditions <- function(loc) {
  
  # model temperature conditions inside a concrete water container (fully shaded,
  # permanent water, 'rock' surface), and compare them to an
  # aboveground, unshaded situation in the same location
  
  # height above the ground (metres) at which to calculate conditions, substrate
  # type, and degree of shade in this place - set to 5cm, rock (ie.
  # concrete) and 100% shade to represent a mosquito resting against side of a
  # concrete water tank
  height_m <- 0.05
  soiltype <- 1
  shade_perc <- 100
  wetness_perc <- 100
  
  micro <- micro_global(
    # place
    loc = loc,
    timeinterval = 365,
    # microclimate characteristics
    Usrhyt = height_m,
    maxshade = shade_perc,
    soiltype = soiltype,
    PCTWET = wetness_perc,
    runmoist = 0
  )
  
  # all outputs, hourly resolution for a year
  list(
    day = seq(0, 365, by = 1/24)[-1],
    habitat = list(
      air_temperature = micro$shadmet[, "TALOC"],
      humidity = micro$shadmet[, "RHLOC"],
      water_temperature = micro$shadsoil[, "D0cm"]
    ),
    outside = list(
      air_temperature = micro$metout[, "TAREF"],
      humidity = micro$metout[, "RH"],
      water_temperature = micro$soil[, "D0cm"]
    )
  )
  
}


# # probability of an adult surviving two weeks at different humidities
# surv_100rh <- function(temp) ds_function(temperature = temp, humidity = rep(100, length(temp))) ^ 14
# surv_50rh <- function(temp) ds_function(temperature = temp, humidity = rep(50, length(temp))) ^ 14
# surv_25rh <- function(temp) ds_function(temperature = temp, humidity = rep(25, length(temp))) ^ 14
# 
# plot(surv_100rh,
#      xlim = c(-20, 100))
# plot(surv_50rh,
#      xlim = c(-20, 100),
#      lty = 2,
#      add = TRUE)
# plot(surv_25rh,
#      xlim = c(-20, 100),
#      lty = 3,
#      add = TRUE)
# min(ds_function(temperature = seq(-20, 100, length.out = 100),
#                  humidity = sample(seq(0, 100, length.out = 100))))

# das_function is a function of temperature and density (larvae per 250ml
# water), the others are functions of temperature

# given our microclimate data, we can now compute these parameters

# For the larval stages (mdr, das) we can use water temperature. The
# experiements use air temperature, but small volumes of water and high humidity
# so that they track with the air temperatures, but in our microclimate there
# will be a buffering effect

conditions <- model_conditions(loc)
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
plot(ds ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "", ylab = "",
     main = "adult survival prob")
plot(das ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     ylim = range(das, das64),
     xlab = "", ylab = "",
     main = "aquatic survival prob \n(low and high density)")
lines(das64 ~ conditions$day,
      lty = 2)
plot(mdr ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "day of year", ylab = "",
     main = "larval development rate")
plot(efd ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "day of year", ylab = "",
     main = "egg laying rate")

# put this into a population dynamic simulation model

# write a function to construct the matrix appropriately, given the larval
# density in the previous timestep

create_matrix <- function(state, water_temperature, mdr, efd, ds, timestep = 1 / 24) {
  
  # given the previous state, and temperature, compute daily das
  das <- das_function(temperature = water_temperature,
                      density = state[1])
  
  # convert all of these to the required timestep (survivals cumuatlive, rates
  # linear)
  das_step <- das ^ timestep
  ds_step <- ds ^ timestep
  mdr_step <- mdr * timestep
  efd_step <- efd * timestep
  
  # construct the matrix
  #     L                A
  # L   das * (1-mdr)    ds * efd
  # A   das * mdr        ds
  matrix(
    c(
      das_step * (1 - mdr_step), # top left
      das_step * (mdr_step), # bottom left
      ds_step * efd_step, # top right
      ds_step # bottom right
    ),
    nrow = 2,
    ncol = 2
  )
  
}

# iterate the state of the model
iterate_state <- function(state, t, water_temperature, mdr, efd, ds) {
  mat <- create_matrix(state = state,
                       water_temperature = water_temperature[t],
                       mdr = mdr[t], efd = efd[t], ds = ds[t])  
  mat %*% state
}

# simulate for a full timeseries, with optional multiple years of burnin
simulate_population <- function(conditions, initial_state = rep(100, 2), burnin_years = 1) {
 
  # add whole year of burnin
  n_times <- length(conditions$water_temperature)
  index <- rep(seq_len(n_times), burnin_years + 1)
  
  # pull out timeseries needed for simulating
  water_temperature <- conditions$water_temperature[index]
  mdr <- mdr_function(conditions$water_temperature[index])
  efd <- efd_function(conditions$air_temperature[index])
  ds <- ds_function(temperature = conditions$air_temperature[index],
                    humidity = conditions$humidity[index])
  
  # simulate the population
  n <- length(index)
  states <- matrix(0, n, 2)
  colnames(states) <- c("aquatic", "adult")
  state <- initial_state
  
  for (t in seq_len(n)) {
    state <- iterate_state(state, t = t,
                           water_temperature = water_temperature,
                           mdr = mdr, efd = efd, ds = ds)
    states[t, ] <- state
  }

  # keep only the final year (post burnin)  
  keep_index <- tail(seq_along(index), n_times)
  states[keep_index, ]  
}


conditions <- model_conditions(loc)
states <- simulate_population(conditions$habitat)

par(mfrow = c(2, 1),
    mar = c(4, 4, 1, 2) + 0.1)
plot(states[, 1] ~ conditions$day,
     ylab = "aquatic stages",
     xlab = "",
     ylim = c(0, max(states[, 1])),
     type = "l")
plot(states[, 2] ~ conditions$day,
     ylab = "adults",
     type = "l",
     ylim = c(0, max(states[, 2])),
     xlab = "day")

summarise_dynamics <- function(states, surface_area_m2 = 10000) {
  
  # given a habitat surace area in square metres, get the population multiplier
  # to scale up the experimental density dependence to get the absolute
  # population sizes for the given pool of water. The experiment used (Evans et
  # al.) has 250ml water in a 'quart size mason jar', which is rather quaint,
  # but not particularly specific. I'm assuming it's a Ball brand 'regular mouth
  # canning jar'. According to masonjarlifestyle.com, that has a 2 3/8" internal
  # diameter. In real money, that's 6.0325cm, for an area of 28.5814687428cm^2,
  # or 0.00285814687m2.
  multiplier <- surface_area_m2 / 0.00285814687
  states <- states * multiplier
  
  # summarise monthly
  hours <- seq(1, 365 * 24)
  timestamps <- as_datetime("2023-01-01 00:00:00") + lubridate::hours(hours)
  months <- lubridate::month(timestamps)
  
  # based on this volume of water, we consider the species cannot persist if
  # there is not at least 1 larva or 1 adult at all times, calculate this for
  # each month (we can summarise annually later and see if all months are
  # suitable)
  min_larvae <- tapply(states[, 1], months, min)
  min_adults <- tapply(states[, 2], months, min)
  
  # calculate the average number of adults present at any given
  # time, in each month
  adult_mean_abundance <- tapply(states[, 2], months, mean)

  data.frame(
    month = 1:12,
    persistence = min_larvae > 0.5 | min_adults > 0.5,
    relative_abundance = adult_mean_abundance
  )  
    
}

# need to summarise with e.g. year-round persistence and average number of
# mosquitoes

# we can calculate adult lifespan and combine with EIP to get an R0 for
# different malarias?

# speed this up? It's fairly fast already

# wrap it up and summarise it at least


# write a wrapper function to run this for a set of pixels, and then profile it
# to make it run faster

calculate_stephensi_suitability <- function(loc) {
  
  # model the microclimates
  conditions <- model_conditions(loc)
  
  results <- list()
  for (microclimate in c("habitat", "outside")) {
      
      # pull out the appropriate conditions and run the simulation
      habitat_conditions <- conditions[[microclimate]]
      population_states <- simulate_population(conditions = habitat_conditions)
    
      # summarise these, and enter into outputs
      summary <- summarise_dynamics(population_states)
      summary$microclimate <- microclimate
      results[[microclimate]] <- summary
      
  }
  
  bind_rows(results$habitat,
            results$outside)
  
}

# almost all of the computation time is in the fortran microclimate model, so
# not worth speeding up
# library(profvis)
# profvis(
#   suit <- calculate_stephensi_suitability(loc)
# )

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
  facet_wrap(~location) +
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
  facet_wrap(~location) +
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
  facet_wrap(~location) +
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

africa_countries <- function () {
  c(
    "AGO",
    "BDI",
    "BEN",
    "BFA",
    "BWA",
    "CAF",
    "CIV",
    "CMR",
    "COD",
    "COG",
    "COM",
    "CPV",
    "DJI",
    "DZA",
    "EGY",
    "ERI",
    "ESH",
    "ETH",
    "GAB",
    "GHA",
    "GIN",
    "GMB",
    "GNB",
    "GNQ",
    "KEN",
    "LBR",
    "LBY",
    "LSO",
    "MAR",
    "MDG",
    "MLI",
    "MOZ",
    "MRT",
    "MUS",
    "MWI",
    "NAM",
    "NER",
    "NGA",
    "RWA",
    "SDN",
    "SEN",
    "SLE",
    "SOM",
    "SSD",
    "STP",
    "SWZ",
    "TCD",
    "TGO",
    "TUN",
    "TZA",
    "UGA",
    "ZAF",
    "ZMB",
    "ZWE"
  )
}

emro_countries <- function () {
  c(
    "AFG",
    "BHR",
    "DJI",
    "EGY",
    "IRN",
    "IRQ",
    "JOR",
    "KWT",
    "LBN",
    "LBY",
    "MAR",
    "PSE",
    "OMN",
    "PAK",
    "QAT",
    "SAU",
    "SOM",
    "SDN",
    "SYR",
    "TUN",
    "ARE",
    "YEM"
  )
}

searo_countries <- function() {
  c(
    "BGD",
    "BTN",
    "PRK",
    "IND",
    "IDN",
    "MDV",
    "MMR",
    "NPL",
    "LKA",
    "THA",
    "TLS"
  )
}

euro_countries_subset <- function() {
  c(
    "ISR"
  )
}

countries <- function() {
  sort(
    unique(
      c(
        africa_countries(),
        emro_countries(),
        searo_countries(),
        euro_countries_subset()
      )
    )
  )
}

countries_exclude <- function() {
  c(
    "PRK",
    "IDN",
    "TLS",
    "THA"
  )
}


region_countries <- function() {
  setdiff(countries(), countries_exclude())
}
# # quick check
# all(region_countries() %in% geodata::country_codes()$ISO3)

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
sfExport(list = list("model_conditions",
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

# need to check model with Mike Kearney, and consder doing an outside water
# source model too

# update the adult survival model to include data from Krajacich et al. (2020)
# https://doi.org/10.1186/s13071-020-04276-y, with an intercept term on study
# and preferentially using that one as the colony is younger (less lab-adapted)
# and survival is longer, more consistent with field observations of long-lived
# An. gambiae s.l. (possibly aestivation)

# Increase the water volume size when determining the persistence suitability,
# because popultions will consist of multiple water tanks


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
  facet_wrap(
    ~location,
    nrow = 1
  ) +
  theme_minimal()


get_rds <- function(url) {
  tf <- tempfile()
  download.file(url, tf)
  readRDS(tf)
}


# read in Whitakker et al extracted data (from forked repo to safeguard against changes)
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
    coords = st_centroid(shape_admin)
  ) %>%
  select(
    -shape_admin1,
    -shape_admin2,
    -shape_admin
  )

whittaker_coords_list <- whittaker_tidied %>%
  pull(coords) %>%
  st_coordinates() %>%
  as_tibble() %>%
  split(seq_len(nrow(.)))

whittaker_results_list <- lapply(whittaker_coords_list, calculate_stephensi_suitability)

index_list <- lapply(whittaker_tidied$id, function(x) tibble(id = x))

whittaker_results <- mapply(bind_cols, whittaker_results_list, index_list, SIMPLIFY = FALSE) %>%
  bind_rows() %>%
  `rownames<-`(NULL)



whittaker_obs_pred <- whittaker_tidied %>%
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
      microclimate == "habitat"
    ),
    by = c("id", month_id = "month")
  ) %>%
  group_by(id) %>%
  mutate(
    abundance = abundance / mean(abundance, na.rm = TRUE),
    relative_abundance = relative_abundance / mean(relative_abundance)
  ) %>%
  # keep only the higher resolution
  filter(!is.na(admin2))


plot_list <- list()
for (country_name in unique(whittaker_obs_pred$country)) {
  
    obs_pred <- whittaker_obs_pred %>%
      filter(country == country_name)
    
    ylims <- range(
      c(obs_pred$abundance,
        obs_pred$relative_abundance)
    )
  
    plot_list[[country_name]] <- obs_pred %>%
      mutate(
        month = factor(month, levels = month.abb)
      ) %>%
      ggplot(
        aes(
          y = abundance,
          x = month,
          group = id
        )
      ) +
      geom_line(
        color = grey(0.6)
      ) +
      geom_point(
        color = grey(0.6)
      ) +
      geom_line(
        aes(
          y = relative_abundance
        )
      ) +
      coord_cartesian(
        ylim = ylims
      ) +
      facet_wrap(~admin2) +
      theme_minimal() +
      ggtitle(country_name)
  
}

# some good fits (in humid areas?), some poor
plot(plot_list$India)

# pretty noisy - lots of variation within admin2s
plot(plot_list$Pakistan)

# peak off by about a month - could be lifespan, as it takes longer to come back
# up from off season?
plot(plot_list$Afghanistan)

# way off
plot(plot_list$Iran)
plot(plot_list$Myanmar)

# data very noisy, one timeseries, not worth using?
plot(plot_list$Djibouti)


