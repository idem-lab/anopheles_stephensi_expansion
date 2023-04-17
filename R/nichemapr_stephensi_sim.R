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

# placename <- "Awash, Ethiopia"
# loc <- c(40.142981, 8.9972474)

# placename <- "Niamey, Niger"
# loc <- c(2.0840132, 13.5127664)

placename <- "Bangui, CAR"
loc <- c(18.561247, 4.365849)

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
