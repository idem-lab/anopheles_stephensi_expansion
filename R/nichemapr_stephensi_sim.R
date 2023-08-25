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
# library(mgcv)

source("R/functions/micro_functions.R")


# load climate-dependent lifehistory functions for An. stephensi
lifehistory_functions_stephensi <- get_lifehistory_functions("An. stephensi")


# pull them out to plot them separately
# daily survival of adults, as a function of temperature and humidity
ds_function <- lifehistory_functions_stephensi$ds_function

# development rate of aquatic stages as a function of water temperature
mdr_function <- lifehistory_functions_stephensi$mdr_function

# daily survival probability of aquatic stages as a function of water
# temperature and density of aquatic stages
das_function <- lifehistory_functions_stephensi$das_function

# daily egg laying as a function of air temperature
efd_function <- lifehistory_functions_stephensi$efd_function

# given our microclimate data, we can now compute these parameters

# For the larval stages (mdr, das) we can use water temperature. The
# experiements use air temperature, but small volumes of water and high humidity
# so that they track with the air temperatures, but in our microclimate there
# will be a buffering effect

addis_loc <- c(38.75960099237313, 8.982497645041477)
awash_loc <- c(40.142981, 8.9972474)
loc <- awash_loc

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
states <- simulate_population(conditions$habitat, lifehistory_functions = lifehistory_functions_stephensi)

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



# comparison against seasonality data

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

whittaker_results_list <- lapply(whittaker_coords_list,
                                 calculate_suitability,
                                 lifehistory_functions = lifehistory_functions_stephensi)

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
  reframe(
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
   


# addis_loc <- c(38.75960099237313, 8.982497645041477)
# loc <- addis_loc
awash_loc <- c(40.142981, 8.9972474)
loc <- awash_loc
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

# to do:

# pull out lifehistory estimation functions into a separate file, and source it

# compute larval density dependence for An. gambiae, from Muriu, or the other
# available data

# make rasters for An. stephensi and An. gambiae, and save monthly densities as
# separate layers

# tidy up visualisation of lifehistory timeseries (ggplot version of modelling
# lifehistory traits and abundances on different timeframes)

# clean up scripts


