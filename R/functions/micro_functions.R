# functions for the microhabitat mecahnistic mdoelling

######
# misc

# download and load an RDS file
get_rds <- function(url) {
  tf <- tempfile()
  download.file(url, tf)
  readRDS(tf)
}

######
# GIS functions (almost certainly replicated elsewhere - so rationalise)

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


######
# lifehistory functions

# reload lifehistory functions from saved objects
rehydrate_lifehistory_function <- function(path_to_object) {
  object <- readRDS(path_to_object)
  do.call(`function`,
          list(object$arguments,
               body(object$dummy_function)))
}

# ensure all values of x are greater than or equal to min
# this is much faster than pmax(x, min)
ensure_gte <- function(x, min = .Machine$double.eps) {
  mask <- as.numeric(x >= min)
  x * mask + min * (1 - mask)
}

# this is much faster than pmax(0, x)
ensure_positive <- function(x) {
  x * as.numeric(x > 0)
}

# this is much faster than pmin(x, max)
enforce_max <- function(x, max) {
  too_big <- x > max
  x[too_big] <- max
  x
}



######
# microclimate modelling

# given a latitude and a longitude, model the temperuture and humidity profile
# at an hourly resolution over an average year, both inside a hypothetical
# concrete water tank, and out in the open air
model_climatic_conditions <- function(loc) {
  
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
    runmoist = 0,
    rainfrac = 0,
    evenrain = 1
  )
  
  # sequence of times, in fractional days of the year
  day <- seq(0, 365, by = 1/24)[-1]
  
  # get smoothed rainfall data
  rainfall <- smooth_rainfall(micro, day)
  
  # all outputs, hourly resolution for a year, with those relevant to either the
  # microclimate or non-microclimate conditions
  list(
    habitat = list(
      day = day,
      air_temperature = micro$shadmet[, "TALOC"],
      humidity = micro$shadmet[, "RHLOC"],
      water_temperature = micro$shadsoil[, "D0cm"],
      windspeed = micro$shadmet[, "VLOC"],
      rainfall = rainfall,
      altitude = micro$elev
    ),
    outside = list(
      day = day,
      air_temperature = micro$metout[, "TAREF"],
      humidity = micro$metout[, "RH"],
      water_temperature = micro$soil[, "D5cm"],
      windspeed = micro$metout[, "VLOC"],
      rainfall = rainfall,
      altitude = micro$elev
    )
  )
  
}

# smooth the nichemapr rainfall output into something more reasonable and smooth
# ona. daily basis
smooth_rainfall <- function(micro, day = NULL) {
  
  # nichemapr gives daily rain per hour, so convert to hourly before smoothing
  rainfall_raw <- micro$RAINFALL / 24
  
  # pull out monthly-constant rainfall and create days variable
  df_fit <- data.frame(
    rainfall = rainfall_raw,
    day = seq(1, 365, length.out = length(rainfall_raw))
  )
  
  # if a new timeseries of days is provided predict to that, otherwise to the
  # training data
  if (is.null(day)) {
    df_pred <- df_fit
  } else {
    df_pred <- data.frame(
      day = day
    )
  }
  
  if (all(rainfall_raw == 0)) {
    return(
      rep(0, length(df_pred$day))
    )
  }
  
  # fit a cyclic spline on log1p transformed data (to handle zero rainfall
  # months)
  m <- mgcv::gam(log1p(rainfall) ~ s(day, bs = "cc", k = 12),
                 data = df_fit,
                 method = "REML")
  
  # predict to the day variable
  pred <- predict(m,
                  newdata = df_pred,
                  type = "link")
  
  # force non-negative
  rainfall_smooth <- pmax(0, expm1(pred))
  
  # ensure totals are correct
  rainfall_smooth <- sum(rainfall_raw) * rainfall_smooth / sum(rainfall_smooth) 
  
  rainfall_smooth
  
}

# extract the modelled monthly rain totals for a location, for plotting
get_rain_months <- function(coords) {
  conditions <- model_climatic_conditions(coords)
  tibble(
    rainfall = conditions$habitat$rainfall,
    date = as.Date("2023-01-01") + conditions$habitat$day
  ) %>%
    mutate(
      # get the month
      month = format(date, "%b"),
    ) %>%
    group_by(
      month
    ) %>%
    summarise(
      rainfall = sum(rainfall),
      .groups = "drop"
    ) %>%
    mutate(
      # get month as a number from 1:12
      month_id = match(month, month.abb),
      # the number for the month after (wrapping around so 13 = 1)
      month_end_id = month_id %% 12 + 1,
      # get the month after
      month_end = month.abb[month_end_id],
      # make both factors
      month = factor(month, levels = month.abb),
      month_end = factor(month_end, levels = month.abb),
    ) %>%
    arrange(month)
}



format_climatic_data <- function(loc) {
  
  # run microclimate model
  climate <- model_climatic_conditions(loc)
  
  # run larval habitat
  ephemeral_habitat_microclimate <- simulate_ephemeral_habitat(
    conditions = climate$habitat)
  ephemeral_habitat_outside <- simulate_ephemeral_habitat(
    conditions = climate$outside)
  
  # format as dataframe
  micro <- climate$habitat %>%
    as_tibble() %>%
    mutate(
      larval_habitat = simulate_ephemeral_habitat(.)
    ) %>% 
    mutate(
      which = "microclimate",
      .before = everything()
    )
  
  outside <- climate$outside %>%
    as_tibble() %>%
    mutate(
      larval_habitat = simulate_ephemeral_habitat(.)
    ) %>% 
    mutate(
      which = "outside",
      .before = everything()
    )
  
  bind_rows(micro, outside) %>%
    select(
      -altitude,
      -windspeed
    ) %>%
    mutate(
      date = as.Date("2023-01-01") + day,
      week = lubridate::week(date),
    ) %>%
    pivot_longer(
      cols = c(
        air_temperature,
        water_temperature,
        humidity,
        rainfall,
        larval_habitat
      ),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(
      variable = case_when(
        variable == "air_temperature" ~ "Air temperature (C)",
        variable == "water_temperature" ~ "Water temperature (C)",
        variable == "humidity" ~ "Humidity (%)",
        variable == "rainfall" ~ "Rainfall (mm)",
        variable == "larval_habitat" ~ "Pooled water\n(relative area)"
      )
    )
  
}

######
# larval habitat modelling

# assume each pool is a cone with 90 degree angle, and define functions to
# convert between measurements
cone_depth_to_radius <- function(depth) {
  depth * sin(pi / 2)
}

cone_depth_to_volume <- function(depth) {
  pi * cone_depth_to_radius(depth) ^ 2 * depth / 3
}

cone_volume_to_depth <- function(volume) {
  ((volume * 3) / (pi * sin(pi / 2) ^ 2)) ^ (1/3)
}

cone_volume_to_surface <- function(volume) {
  # depth <- cone_volume_to_depth(volume)
  # pi * cone_depth_to_radius(depth) ^ 2
  # with the right angle asusmption, we have sin(pi / 2) = 1, which cancels
  # leaving this (which is a bit faster to evaluate in the larval habitat model)
  # pi * (3 / pi) ^ (2/3) * volume ^ (2/3)
  3.046474 * volume ^ 0.666667
}

# use the Buck equation to estimate the vapour pressure in kPa
# see: https://en.wikipedia.org/wiki/Vapour_pressure_of_water
# and: https://en.wikipedia.org/wiki/Arden_Buck_equation
saturation_vapour_pressure_kpa <- function(temperature_celsius) {
  0.61121 * exp((18.678 - temperature_celsius / 234.5) * (temperature_celsius / (257.14 + temperature_celsius)))
}

# air pressure as a function of altitude
air_pressure_kpa <- function(altitude_metres) {
  # At low altitudes above sea level, the pressure decreases by about 1.2 kPa
  # (12 hPa) for every 100  metres.
  p0 <- 101.325 # kPa at sea level
  p0 - 1.2 * altitude_metres / 100
}

# compute the amount of evaporation from a water body (kg per hour) as a function of
# temperature, humidity, altitude, windspeed, and surface area
evaporation_rate_kg_h <- function(temperature_c,
                                  relative_humidity,
                                  altitude_m,
                                  windspeed_m_s,
                                  surface_area_m2) {
  
  # https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html
  # rate of evaporation gh (kg/hour) from a water body is:
  # gh = theta * A * (xs - x)
  # theta = (25 + 19 v) = evaporation coefficient (kg/m2h)
  # v = velocity of air above the water surface (m/s)
  # A = surface_area (m^2)
  # xs = *maximum* humidity ratio of saturated air at the same temperature as
  #   the water surface (kg/kg)  (kg H2O in kg Dry Air)
  # x = humidity ratio of the air
  
  v <- windspeed_m_s
  theta <- (25 + 19 * v)
  A <- surface_area_m2
  
  # compute the maximum humidity ratio xs, given the saturation vapour pressure p_ws
  # and atmospheric pressure p_a
  # xs = 0.62198 p_ws / (p_a - p_ws)
  p_ws <- saturation_vapour_pressure_kpa(temperature_c)
  p_a <- air_pressure_kpa(altitude_m)
  xs <- 0.62198 * p_ws / (p_a - p_ws)
  
  # now compute the humidity ratio x for the actual air, from the air pressure and
  # actual vapour pressure
  
  # relative humidity is the ratio of the partial pressure of water vapour in the
  # air (p), to the saturation vapour pressure (p_s), so use it to convert back to
  # the partial pressure of water vapour p
  # RH = 100 * p_w / p_ws
  p_w <- p_ws * (relative_humidity / 100)
  x <- 0.62198 * p_w / (p_a - p_w)
  
  # put it all together
  evaporation_kg_hour <- theta * A * (xs - x)
  evaporation_kg_hour
}

# iterate the water volume in the cone, with accounting for evaporation and
# rainfall
iterate_cone_volume <- function(water_volume,
                                t,
                                rainfall_mm_h,
                                evaporation_rate_kg_h_m2,
                                # temperature_c,
                                # relative_humidity,
                                # altitude_m,
                                # windspeed_m_s,
                                max_cone_volume = pi/3,
                                catchment_area = pi) {
  
  # at each 1h timestep, compute evaporation before inflow
  surface_area_m2 <- cone_volume_to_surface(water_volume)
  loss_kg_h <- evaporation_rate_kg_h_m2[t] * surface_area_m2
  
  # this is already a 1h timestep, so just convert this to m^3
  # 1kg = 1,000cm3
  # 1m3 = 1,000,000 cm3, so
  # 1m3 = 1,000kg
  loss <- loss_kg_h / 1000
  
  # lose water, avoiding negative volume
  water_volume <- ensure_positive(water_volume - loss)
  
  # 1000mm per metre, and catchment area is in m2, so convert to m3 (including
  # multiplier on inflow)
  gain_m_h <- rainfall_mm_h[t] / 1000
  gain <- catchment_area * gain_m_h
  
  # gain water, capping maximum volume
  water_volume <- enforce_max(water_volume + gain, max = max_cone_volume)
  
  water_volume
  
}

# given timeseries of conditions (including rainfall, windspeed and altitude),
# compute a timeseries of water surface area from a simple cone model
simulate_ephemeral_habitat <- function(conditions,
                                       initial_volume = 0,
                                       burnin_years = 1,
                                       max_cone_depth = 1,
                                       inflow_multiplier = 1) {
  
  # add whole year of burnin
  n_times <- length(conditions$water_temperature)
  index <- rep(seq_len(n_times), burnin_years + 1)
  
  # pull out timeseries needed for simulating
  rainfall <- conditions$rainfall[index]
  air_temperature <- conditions$air_temperature[index]
  relative_humidity <- conditions$humidity[index]
  windspeed <- conditions$windspeed[index]

  # and fixed info  
  altitude <- conditions$altitude

  # set up the cone model
  
  # define the cone - assume it cone has a maximum depth of 1m, so a maximum
  # volume of 1.05m3, beyond which it cannot get more full
  max_cone_volume <- cone_depth_to_volume(max_cone_depth)
  
  # the catchment area is arbitrary (modelling only relative
  # abundance), but set it to the maximum cone surface area - pi!
  catchment_area <- inflow_multiplier * cone_volume_to_surface(max_cone_volume)
  
  
  # precalculate the evaporation rate per unit surface area
  loss_kg_h_m2 <- evaporation_rate_kg_h(
    temperature_c = air_temperature,
    relative_humidity = relative_humidity,
    altitude_m = altitude,
    windspeed_m_s = windspeed,
    surface_area_m2 = 1
  )
    
  # simulate the water volume
  n <- length(index)
  volumes <- numeric(length = length(index))
  volume <- 0
  
  for (t in seq_along(volumes)) {
    volume <- iterate_cone_volume(water_volume = volume,
                                  t = t,
                                  rainfall_mm_h = rainfall,
                                  evaporation_rate_kg_h_m2 = loss_kg_h_m2,
                                  # temperature_c = air_temperature,
                                  # relative_humidity = relative_humidity,
                                  # altitude_m = altitude,
                                  # windspeed_m_s = windspeed,
                                  max_cone_volume = max_cone_volume,
                                  catchment_area = catchment_area)
    volumes[t] <- volume
  }

  # keep only the final year (post burnin)  
  keep_index <- tail(seq_along(index), n_times)
  volumes <- volumes[keep_index]  
  
  # convert to (depth then) surface area
  surface_areas <- cone_volume_to_surface(volumes)
  
}




######
# population modelling

# Build a simple two stage model (aquatic stages and adults), with effects of
# density dependence (daily aquatic survival; DAS), water temperature (DAS and
# aquatic development rate; MDR), air temperature (adult survival; DS, and egg
# laying; EFD), and humidity (DS). Construct as a dynamic matrix.

# construct the matrix appropriately, given the larval
# density in the previous timestep
create_matrix <- function(state,
                          das_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          ds,
                          timestep = 1 / 24) {
  
  # given the previous state, surface area, and temperature, compute daily
  # aquatic survival (density- and temperature-dependent)
  das <- das_function(temperature = water_temperature,
                      density = state[1] / larval_habitat_area)
  
  # convert all of these to the required timestep (survivals cumulative, rates
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
iterate_state <- function(state,
                          t,
                          das_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          ds) {
  mat <- create_matrix(state = state,
                       das_function,
                       larval_habitat_area = larval_habitat_area[t],
                       water_temperature = water_temperature[t],
                       mdr = mdr[t],
                       efd = efd[t],
                       ds = ds[t])  
  mat %*% state
}

# access a list of the lifehistory functions needed for the named species
get_lifehistory_functions <- function(
    species = c("An. stephensi", "An. gambiae"),
    storage_path = "data/life_history_params/dehydrated") {
  
  # enforce the species label
  species <- match.arg(species)
  
  # load the daily adult survival for either An. gambiae or An. stephensi
  ds_temp_humid = rehydrate_lifehistory_function(
    file.path(storage_path, "ds_temp_humid.RDS")
  ) 
  
  # subset to this species
  ds_function <- function(temperature, humidity) {
    ds_temp_humid(temperature, humidity, species = species)
  }
  
  # for the others, load the relevant RDS object
  species_suffix <- switch(species,
                           "An. stephensi" = "As",
                           "An. gambiae" = "Ag")
  
  # development rate of aquatic stages as a function of water temperature
  mdr_function = rehydrate_lifehistory_function(
    file.path(storage_path,
              sprintf("mdr_temp_%s.RDS", species_suffix))
  )
  
  # daily survival probability of aquatic stages as a function of water
  # temperature and density of aquatic stages
  das_function = rehydrate_lifehistory_function(
    file.path(storage_path,
              sprintf("das_temp_dens_%s.RDS", species_suffix))
  )
  
  # daily egg laying as a function of air temperature
  efd_function = rehydrate_lifehistory_function(
    file.path(storage_path,
              sprintf("efd_temp_%s.RDS", species_suffix))
  )
  
  # return as a named list of functions
  list(
    ds_function = ds_function,
    mdr_function = mdr_function,
    das_function = das_function,
    efd_function = efd_function
  )
  
}

# simulate for a full timeseries, with optional multiple years of burnin
simulate_population <- function(
    conditions,
    lifehistory_functions,
    larval_habitat_area = rep(pi, length(conditions$day)),
    initial_state = rep(100, 2),
    burnin_years = 1) {
  
  # add whole year of burnin
  n_times <- length(conditions$water_temperature)
  index <- rep(seq_len(n_times), burnin_years + 1)
  
  # pull out timeseries needed for simulating
  water_temperature <- conditions$water_temperature[index]
  mdr <- lifehistory_functions$mdr_function(
    conditions$water_temperature[index])
  efd <- lifehistory_functions$efd_function(
    conditions$air_temperature[index])
  ds <- lifehistory_functions$ds_function(
    temperature = conditions$air_temperature[index],
    humidity = conditions$humidity[index])
  larval_habitat_area <- larval_habitat_area[index]
  
  # simulate the population
  n <- length(index)
  states <- matrix(0, n, 2)
  colnames(states) <- c("aquatic", "adult")
  state <- initial_state
  
  # pass int he daily aquatic survival function, as it is the only one that is
  # dynamic (depends on the previous state, so introduces density dependence)
  for (t in seq_len(n)) {
    state <- iterate_state(state,
                           t = t,
                           das_function = lifehistory_functions$das_function,
                           larval_habitat_area = larval_habitat_area,
                           water_temperature = water_temperature,
                           mdr = mdr,
                           efd = efd,
                           ds = ds)
    states[t, ] <- state
  }
  
  # keep only the final year (post burnin)  
  keep_index <- tail(seq_along(index), n_times)
  states[keep_index, ]  
}

summarise_dynamics <- function(states, surface_area_m2 = 10000) {
  
  # given a habitat surface area in square metres, get the population multiplier
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

# given a single set of coordinates and a set of lifehistory functions (relating
# microclimatic conditions to lifehistory parameters), return monthly
# suitability measures for that species in all combinations of the
# (micro)climate and larval habitat types
calculate_suitability <- function(loc,
                                  lifehistory_functions,
                                  microclimates = c("habitat", "outside"),
                                  larval_habitats = c("permanent", "ephemeral")) {
  
  # check the microclimates and larval habitats to model are valid options
  if(!all(microclimates %in% c("habitat", "outside"))) {
    stop("invalid microclimates specified")
  }
  
  if(!all(larval_habitats %in% c("permanent", "ephemeral"))) {
    stop("invalid larval habitats specified")
  }
  
  # model the microclimates
  climatic_conditions <- model_climatic_conditions(loc)
  
  results <- list()
  for (microclimate in microclimates) {
    
    # pull out the appropriate microclimatic conditions
    microclimate_conditions <- climatic_conditions[[microclimate]]
    
    for (larval_habitat in larval_habitats) {
    
      # compute the appropriate area of larval habitat  
      if (larval_habitat == "ephemeral") {
        
        # simulate ephemeral larval habitat under these local climatic conditions
        larval_habitat_area <- simulate_ephemeral_habitat(
          conditions = microclimate_conditions)

      } else {
        # use pi for a permanent water body because it's the survace area (in m2) of a full, 1m deep
        # right-angle cone. But the value doesn't matter because we are only
        # after relative abundances
        larval_habitat_area <- rep(pi, length(microclimate_conditions$day))
      }
      
      # and run the simulation
      population_states <- simulate_population(
        conditions = microclimate_conditions,
        lifehistory_functions = lifehistory_functions,
        larval_habitat_area = larval_habitat_area)
      
      # summarise these, and enter into outputs
      summary <- summarise_dynamics(population_states)
      summary$microclimate <- microclimate
      summary$larval_habitat <- larval_habitat
      
      results[[microclimate]][[larval_habitat]] <- summary
      
    }
    
  }

  # combine them all and return  
  bind_rows(results$habitat$permanent,
            results$habitat$ephemeral,
            results$outside$permanent,
            results$outside$ephemeral)
  
}



