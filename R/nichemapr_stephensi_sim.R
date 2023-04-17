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
  PCTWET = wetness_perc
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



# need to plug these into the temperature suitability model


# To do:

# find all the centroids of the climate raster, for the Afro and EMRO regions

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

# port the microclimate conditions into the adult survival model
# load the mgcv survival model 
library(mgcv)
survival_model <- readRDS("data/survival_model/survival_model.RDS")

survival <- function (temperature = 20,
                      humidity = 100,
                      model = NULL,
                      use_wild_correction = FALSE,
                      period_days = 1) {
  
  if (is.null(model)) {
    stop ("model not specified")
  }
  
  # construct dataframe for prediction
  df_pred <- data.frame(temperature = temperature,
                        humidity = humidity,
                        sex = "F",
                        off = 0)
  
  # predict cumulative survival probability after one day
  survival <- 1 - predict(model, df_pred, type = "response")
  
  if (use_wild_correction) {
    wild_correction <- attr(model, "wild_correction")
    survival <- survival * wild_correction
  }
  
  # correct to given time period over which to calculate survival
  survival <- survival ^ period_days
  
  survival
  
}


# get hourly survival probabilities
survival_microclimate_all <- survival(temperature = microclimate_temperature_all,
                                      humidity = microclimate_humidity_all,
                                      model = survival_model,
                                      use_wild_correction = TRUE,
                                      period_days = 1 / 24)

survival_outside_all <- survival(temperature = outside_temperature_all,
                                 humidity = outside_humidity_all,
                                 model = survival_model,
                                 use_wild_correction = TRUE,
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

# plug the survival curves (and other An. stephensi temperature parameters) into
# the temperature suitability function(s)

# model the duration of the gonotrophic cycle for An stephensi, from Shapiro et
# al. using the enzyme kinetics model of Focks et al.

paaijmans_mean_gono <- tibble::tribble(
  ~temp, ~days,
  22, 5.2,
  24, 5.0,
  26, 4.1
) %>%
  mutate(
    study = "Paaijmans"
  )
shapiro_mean_gono <- readxl::read_excel(
  "data/ovi_data/PLoS.Biology.DataFigures.Supplemental.xlsx",
  sheet = "Fig.4",
  range = "A1:B7") %>%
  rename(
    temp = Temp,
    days = `G-cycle`
  ) %>%
  mutate(
    temp = str_remove(temp, pattern = "ºC"),
    temp = as.numeric(temp),
    study = "shapiro"
  )

mean_gono <- bind_rows(
  paaijmans_mean_gono,
  shapiro_mean_gono
)

# fit a Bayesian model of degree day accumulation for the first gonotrophic
# cycle
library(greta)

# based on eq from Brady et al., gives the development rate per hour at the
# given temperature
focks_gono_rate <- function(t_c, rho, h_a, h_h, t_0_5) {
  # convert to kelvin
  t_k <- t_c + 273.15
  # compute equation components
  a <- rho * t_k / 298
  b <- exp((h_a / 1.987) * ((1 / 298) - (1 / t_k)) )
  c <- exp((h_h / 1.987) * ((1 / t_0_5) - (1 / t_k)) )
  # return gonotrophic cycle development rate per hour
  a * b / (1 - c)
}

# take the prior means and initial values from Sharpe & de Michelle (via Focks)
# estimates for Ae. aegypti
rho_mean <- 0.00898
h_a_kcal_mean <- 15725.23 / 1000
h_h_mcal_mean <- 1756481.07 / 1000000
t_0_5_mean <- 447.17

rho <- normal(rho_mean, 1, truncation = c(0, Inf))
h_a_kcal <- normal(h_a_kcal_mean, 100, truncation = c(0, Inf))
h_h_mcal <- normal(h_h_mcal_mean, 10, truncation = c(0, Inf))
t_0_5 <- normal(t_0_5_mean, 10, truncation = c(0, Inf))
obs_sd <- normal(0, 1, truncation = c(0, Inf))

h_a <- h_a_kcal * 1000
h_h <- h_h_mcal * 1000000

dev_rates <- focks_gono_rate(shapiro_mean_gono$temp,
                             rho = rho,
                             h_a = h_a,
                             h_h = h_h,
                             t_0_5 = t_0_5)

# convert these to the expected number of days
expected_days <- 1 / (dev_rates * 24)
distribution(shapiro_mean_gono$days) <- normal(expected_days, obs_sd)
m <- model(rho, h_a, h_h, t_0_5, obs_sd)

n_chains <- 4
inits <- replicate(n_chains,
             initials(
               rho = rho_mean,
               h_a_kcal = h_a_kcal_mean,
               h_h_mcal = h_h_mcal_mean,
               t_0_5 = t_0_5_mean
             ),
             simplify = FALSE)

draws <- mcmc(m, chains = n_chains,
              initial_values = inits)

# pull out the parameter values to use as priors
plot(draws)
summary(draws)

range(shapiro_mean_gono$temp)
new_temps <- seq(-20, 60, length.out = 100)
new_gc_rates <- focks_gono_rate(new_temps,
                           rho = rho,
                           h_a = h_a,
                           h_h = h_h,
                           t_0_5 = t_0_5)

new_gcs <- 1 / (new_gc_rates * 24)

post_means <- as.list(summary(draws)$statistics[, "Mean"])
sims <- calculate(new_gcs,
                  new_gc_rates,
                  values = draws, nsim = 1000)

preds <- sims$new_gcs[, , 1]
pred_means <- colMeans(preds)
pred_quants <- apply(preds, 2, quantile, c(0.025, 0.975))
par(mfrow = c(1, 1),
    mar = c(5, 5, 2, 2))
plot(pred_means ~ new_temps,
     type = "n",
     ylim = c(0, 9),
     xlim = c(18, 38),
     axes = FALSE,
     xlab = "Temperature (celsius)",
     ylab = "Mean gonotrophic cycle (days)")
polygon(x = c(new_temps, rev(new_temps)),
       y = c(pred_quants[1, ], rev(pred_quants[2, ])),
       lty = 0,
       col = grey(0.9))
lines(pred_means ~ new_temps, pch = 16)
axis(1)
axis(2, las = 2)
points(shapiro_mean_gono$days ~ shapiro_mean_gono$temp,
       col = "red", pch = 15)
box()

# create a spline interpolation of this function for use in the cohort
# simulation
pred_rates <- colMeans(sims$new_gc_rates[, , 1])
gc_function <- splinefun(x = new_temps,
                         y = pred_rates)

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


# install and load Oli's tempsuit package, and wrap it up in a more sane
# interface so it does't constantly segfault and overwrite objects

# install.packages("data/tempsuitcalc_0.1.tar.gz", type = "source", repos = NULL)

library(tempsuitcalc)

# create survival
myP_vec <- survival(temperature = temperature,
                    humidity = humidity,
                    model = survival_model, 
                    use_wild_correction = TRUE,
                    period_days = 1 / 12)

# replicate into matrix (we don't have an age effect in this survival model)
myP <- myP_vec %*% matrix(1, 1, max_length)

# write a function to more safely handle the inputs and outputs
simulate_suitability <- function(temperature_hourly,
                            humidity_hourly,
                            pathogen = c("vivax", "falciparum")) {

  pathogen <- match.arg(pathogen)
  
  # create indices to pull out 2h bin values, with burnin and burnout  
  bin_index_2h_year <- seq(1, 365 * 24, by = 2)
  bin_index_2h_burnout <- seq(1, 5 * 7 * 24 - 24, by = 2)
  bin_index_2h_burnin <- 365 * 24 - rev(seq(1, 5 * 7 * 24, by = 2)) - 2
  bin_index_2h <- c(bin_index_2h_burnin, bin_index_2h_year, bin_index_2h_burnout)
  max_length <- as.integer(10416 / 2)
  stopifnot(identical(max_length, length(bin_index_2h)))
  
  # pull out full timeseries
  temperature <- temperature_hourly[bin_index_2h]
  humidity <- humidity_hourly[bin_index_2h]
  
  # compute life history parameters over time
  myDD <- eip_function(temperature, pathogen = pathogen)
  myovi <- gc_function(temperature)
  
  # create survival
  myP_vec <- survival(temperature = temperature,
                      humidity = humidity,
                      model = survival_model, 
                      use_wild_correction = TRUE,
                      period_days = 1 / 12)
  
  # replicate into matrix (we don't have an age effect in this survival model)
  myP <- myP_vec %*% matrix(1, 1, max_length)
  
  # create dummy output vectors
  Zout <- rep(0, max_length)
  Zoutovi <- rep(0, max_length)
  
  # run simulation (overwrites Zout and Zoutovi in place :0 )
  tempsuitcalc::cohort.simulate(myDD = myDD,
                                myovi = myovi,
                                myP = myP,
                                Zout = Zout,
                                Zoutovi = Zoutovi)
  
  # get index to middle years, and return only these
  index <- seq_len(365 * 12) + 5 * 7 * 12
  
  list(Zout = Zout[index],
       Zoutovi = Zoutovi[index])
  
}

summarise_suitability <- function(Zout){
  
  # summarise suitability (number/presence of infectious/reproductive
  # individuals)
  Zout_total <- sum(Zout)
  NID <- (Zout > 0) / 12
  NID_total <- sum(NID)
  
  # sum over weeks, and assign excess days into the last week
  weeks <- c(rep(seq_len(52), each = 7*12), rep(52, (365 - 52 * 7) * 12))
  Zout_weekly <- tapply(Zout, weeks, FUN = sum)
  NID_weekly <- tapply(NID, weeks, FUN = sum)
  
  return(
    list(
      Zout_total = Zout_total,
      NID_total = NID_total,
      Zout_weekly = Zout_weekly,
      NID_weekly = NID_weekly
    )
  )
}

monthly_means <- function(Zout) {
  # aggregate outputs by month
  hours <- seq(1, 365 * 24, by = 2)
  timestamps <- as_datetime("2023-01-01 00:00:00") + lubridate::hours(hours)
  months <- lubridate::month(timestamps)
  tapply(Zout, months, FUN = mean)
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
    PCTWET = wetness_perc
  )
  
  # all outputs, hourly resolution for a year
  list(
    habitat = list(
      temperature = micro$shadmet[, "TALOC"],
      humidity = micro$shadmet[, "RHLOC"]
    ),
    outside = list(
      temperature = micro$metout[, "TAREF"],
      humidity = micro$metout[, "RH"]
    )
  )
  
}

# write a wrapper function to run this for a set of pixels, and then profile it
# to make it run faster

calculate_stephensi_suitability <- function(loc) {
  
  # model the microclimates
  conditions <- model_conditions(loc)
  
  results <- list()
  for (microclimate in c("habitat", "outside")) {
      
      habitat_conditions <- conditions[[microclimate]]
      
      sim_vivax <- simulate_suitability(
        temperature_hourly = habitat_conditions$temperature,
        humidity_hourly = habitat_conditions$humidity,
        pathogen = "vivax" 
      )
      
      sim_falciparum <- simulate_suitability(
        temperature_hourly = habitat_conditions$temperature,
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


# work out why it's so optimistic about outside the habitat - is it because reduced survival
# doesn't have as much impact as the speed of the gonotrophic cycle or EIP?

# different survival is definitely being passed through.

# maybe the repeated releases are nullifying survival?





# run this over all pixels to map year-round microclimate (and outside
# microclimate) suitability







