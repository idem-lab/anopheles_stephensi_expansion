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

placename <- "Awash, Ethiopia"
loc <- c(40.142981, 8.9972474)

# placename <- "Niamey, Niger"
# loc <- c(2.0840132, 13.5127664)

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

# load the mgcv survival model 
library(mgcv)
survival_model <- readRDS("data/survival_model/survival_model.RDS")

# function to compute adult survival at temperature and humidity combinations
survival <- function (temperature = 20,
                      humidity = 100,
                      use_wild_correction = TRUE,
                      period_days = 1) {
  
  # construct dataframe for prediction
  df_pred <- data.frame(temperature = temperature,
                        humidity = humidity,
                        sex = "F",
                        off = 0)
  
  # predict cumulative survival probability after one day
  survival <- 1 - predict(survival_model, df_pred, type = "response")
  
  if (use_wild_correction) {
    wild_correction <- attr(survival_model, "wild_correction")
    survival <- survival * wild_correction
  }
  
  # correct to given time period over which to calculate survival
  survival <- survival ^ period_days
  
  survival
  
}


# get hourly survival probabilities
survival_microclimate_all <- survival(temperature = microclimate_temperature_all,
                                      humidity = microclimate_humidity_all,
                                      period_days = 1 / 24)

survival_outside_all <- survival(temperature = outside_temperature_all,
                                 humidity = outside_humidity_all,
                                 period_days = 1 / 24)

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
      col = "dark green", lwd = 2)

# function for EIP of pathogens
eip_function <- function(temp, pathogen = c("falciparum", "vivax"), period_days = 1 / 12) {
  pathogen <- match.arg(pathogen)
  # use degree day model as in Gething:
  # https://doi.org/10.1186/1756-3305-4-92
  phi <- switch(pathogen,
                vivax = 105,
                falciparum = 111)
  temp_min <- switch(pathogen,
                  vivax = 14.5,
                  falciparum = 16)
  
  degree_days <- pmax(0, temp - temp_min) * period_days

  # scale it 1 is full development
  degree_days / phi
}


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
      humidity = micro$metout[, "RH"]
    )
  )
  
}

# the temperature suitability model is no good in this context, because dumping
# 100 mosquitoes per day into the environment makes the EIP effect more dominant
# than the survival effect. Ie. it doesn't tell us about population persistence.
# We need actual population dynamics.

# Build a simple two stage model (aquatic stages and adults), with aquatic stage
# density dependence and possibly temperature (aquatic survival), adult
# temperature and humidity survival (adult survival) larval development (aquatic
# -> adult rate) as a function of water temperature, and gonotrophic cycle
# (adult -> aquatic rate) as a function of air temperature. Construct as a
# dynamic matrix.

# we have the adult humidity and survival combinations, and we can get the other
# parameters from a combination of Villena et al.
# https://doi.org/10.1002/ecy.3685 and Evans et al.

# we'll need to reestimate the curves from the datasets, as the posteriors and
# parameter estimates are not provided for stephensi *facepalm emoji*

# for the rate of development from eggs to larvae, we want the mosquito
# development rate 'MDR' from Villena et al. Villena (and Johnson 2015) don't
# define MDR, but cite Mordecai 2013 which defines it as 'larval development
# rate' and refers to Bayoh and Lindsay. Bayoh and lindsay do a number of
# different things, but it's clear from the data in Villena that they are using
# the time from egg to adult, so MDR is the inverse of the full duration of the
# aquatic stage, which is what we want.

# We need to re-estimate parameters for MDR (larval development, or rate of
# moving from the larval to adult stages, in days; asymmetric), PEA (proportion
# of eggs surviving to adulthood; concave down), and EFD (eggs per female per
# day; concave down). We then use the data from Evans to compute the rate of
# decline of PEA (on the logit scale) as a function of larval density. These are
# all computed in "R/refit_villena_lifehistory.R"

mdr_function_data <- readRDS("data/life_history_params/mdr_function_data.RDS")
pea_function_data <- readRDS("data/life_history_params/pea_function_data.RDS")
efd_function_data <- readRDS("data/life_history_params/efd_function_data.RDS")

mdr_function <- function(temperature) {
  raw_fun <- with(mdr_function_data, splinefun(temperature, value))
  pmax(0, raw_fun(temperature))
}  

efd_function <- function(temperature) {
  raw_fun <- with(efd_function_data, splinefun(temperature, value))
  pmax(0, raw_fun(temperature))
}  

# need to do 2d spline function
pea_function <- function(temperature, density) {
  raw_fun <- with(pea_function_data, splinefun(temperature, raw_value))
  raw_vals <- pmax(0, raw_fun(temperature))
  plogis(qlogis(raw_vals) + density * pea_function_data$density_coef[1])
}  

# plot these and check they don't produce negatives
plot(mdr_function, xlim = c(-20, 100))
min(mdr_function(seq(-20, 100, by = 0.01)))

plot(efd_function, xlim = c(-20, 100))
min(efd_function(seq(-20, 100, by = 0.01)))

pea_0d <- function(temp) pea_function(temp, density = 0)
pea_32d <- function(temp) pea_function(temp, density = 32)
pea_64d <- function(temp) pea_function(temp, density = 64)
plot(pea_0d,
     xlim = c(-20, 100))
plot(pea_32d,
     xlim = c(-20, 100),
     lty = 2,
     add = TRUE)
plot(pea_64d,
     xlim = c(-20, 100),
     lty = 3,
     add = TRUE)
min(pea_function(temperature = seq(-20, 100, length.out = 100),
                 density = sample(seq(0, 100, length.out = 100))))

# probability of surviving two weeks at different humidities
surv_100rh <- function(temp) survival(temperature = temp, humidity = rep(100, length(temp))) ^ 14
surv_50rh <- function(temp) survival(temperature = temp, humidity = rep(50, length(temp))) ^ 14
surv_25rh <- function(temp) survival(temperature = temp, humidity = rep(25, length(temp))) ^ 14

plot(surv_100rh,
     xlim = c(-20, 100))
plot(surv_50rh,
     xlim = c(-20, 100),
     lty = 2,
     add = TRUE)
plot(surv_25rh,
     xlim = c(-20, 100),
     lty = 3,
     add = TRUE)
min(survival(temperature = seq(-20, 100, length.out = 100),
                 humidity = sample(seq(0, 100, length.out = 100))))

# how to scale the proportion of eggs surviving to adulthood, as a function of
# temperature, since we apply mortality on each (hourly) timestep - use the
# expected average duration of the stage at that temperature

pea_daily_function <- function(temperature, density) {
  pea_overall <- pea_function(temperature, density)
  expected_duration_days <- 1 / mdr_function(temperature)
  pea_daily <- pea_overall ^ mdr_function(temperature)
  pea_daily
}

# pea_function is a function of temperature and density (larvae per 250ml
# water), the others are functions of temperature

# given our microclimate data, we can now compute these parameters

# For the larval stages (mdr, pea) we can use water temperature. The
# experiements use air temperature, but small volumes of water and high humidity
# so that they track with the air temperatures, but in our microclimate there
# will be a buffering effect

conditions <- model_conditions(loc)
scaling <- 1/24
mdr <- mdr_function(conditions$habitat$water_temperature) ^ scaling
pea <- pea_daily_function(conditions$habitat$water_temperature, density = 0) ^ scaling
pea64 <- pea_daily_function(conditions$habitat$water_temperature, density = 64) ^ scaling
efd <- efd_function(conditions$habitat$air_temperature) ^ scaling
surv <- survival(temperature = conditions$habitat$air_temperature,
                 humidity = conditions$habitat$humidity) ^ scaling




par(mfrow = c(2, 2),
    oma = rep(0, 4),
    mar = c(3, 3, 2, 1) + 0.1)
plot(surv ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "", ylab = "",
     main = "adult survival (hourly)")
plot(pea ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "", ylab = "",
     main = "aquatic survival (hourly)")
lines(pea64 ~ conditions$day,
      lty = 2)
plot(mdr ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "day of year", ylab = "",
     main = "larval development rate (hourly)")
plot(efd ~ conditions$day,
     type = "l",
     xlim = c(10, 14),
     xlab = "day of year", ylab = "",
     main = "egg laying (hourly)")

# why is aquatic survival so janky?
plot(pea_function(conditions$habitat$water_temperature,
                  density = 0) ~ conditions$day,
     xlim = c(10, 14), type = "l")

library(tidyverse)

plot(conditions$habitat$water_temperature ~conditions$day,
     xlim = c(10, 14), type = "l")


# this is killing larvae when the temperature temporarily drops below T0, how to
# avoid this? degree day type model modification? Reparameterise without using
# hard boundaries?





# write a wrapper function to run this for a set of pixels, and then profile it
# to make it run faster

calculate_stephensi_suitability <- function(loc) {
  
  # model the microclimates
  conditions <- model_conditions(loc)
  
  results <- list()
  for (microclimate in c("habitat", "outside")) {
      
      habitat_conditions <- conditions[[microclimate]]
      
      sim_vivax <- simulate_suitability(
        temperature_hourly = habitat_conditions$air_temperature,
        humidity_hourly = habitat_conditions$humidity,
        pathogen = "vivax" 
      )
      
      sim_falciparum <- simulate_suitability(
        temperature_hourly = habitat_conditions$air_temperature,
        humidity_hourly = habitat_conditions$humidity,
        pathogen = "falciparum" 
      )
      
      # summarise these, and enter into outputs
      
      results[[microclimate]] <- tibble(
        microclimate = microclimate,
        month = format(ISOdate(2023,1:12,1),"%B"),
        vivax_suitability = monthly_means(sim_vivax$Zout),
        falciparum_suitability = monthly_means(sim_falciparum$Zout),
        reproduction_suitability = monthly_means(sim_vivax$Zoutovi),
        vivax_NID = monthly_means(sim_vivax$Zout > 0),
        falciparum_NID = monthly_means(sim_falciparum$Zout > 0),
        reproduction_NID = monthly_means(sim_vivax$Zout > 0)
      )
      
  }
  
  bind_rows(results$habitat,
            results$outside)
  
}
  
suit <- calculate_stephensi_suitability(loc)

library(ggplot2)

suit %>%
  mutate(
    month = factor(month, levels = unique(month))
  ) %>%
  ggplot(
    aes(x = month,
        y = falciparum_NID,
        group = microclimate,
        colour = microclimate)
  ) +
  geom_line()




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


# aggregate it 9x to get a lower resolution set of coordinates for first pass -
# original resolution is ~18.5km at the equator
region_raster_mask_agg < aggregate(region_raster_mask, 9)
coords_agg <- xyFromCell(region_raster_mask_agg, cells(region_raster_mask_agg))







