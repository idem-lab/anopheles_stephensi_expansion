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

# run this over all pixels to map year-round microclimate (and outside
# microclimate) suitability







