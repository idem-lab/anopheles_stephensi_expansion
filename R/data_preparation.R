library(dplyr)
library(readxl)
library(lubridate)
library(readr)

# data from https://apps.who.int/malaria/maps/threats
mtm_raw <- read_xlsx(
  path = "data/MTM_INVASIVE_VECTOR_SPECIES_20230303.xlsx",
  sheet = "Data"
)


pa_data <- mtm_raw %>%
  mutate(
    #date = paste(YEAR_END, MONTH_END) %>% ym
    # date = case_when(
    #   !is.na(MONTH_START) ~ paste(YEAR_START, MONTH_START, "01") %>% ymd,
    #   !is.na(YEAR_START) ~  paste(YEAR_START, "JAN", "01") %>% ymd,
    #   TRUE ~ NA_Date_
    # ), # something sus going on
    year = YEAR_START,
    presence = if_else(INVASIVE_STATUS == "not found", 0, 1), # either listed as "invasive", or "native" if present, or "not found" if absent
    presence = if_else(is.na(presence), 1, presence)
  ) %>%
  select(
    #date,
    year,
    presence,
    x = LONGITUDE,
    y = LATITUDE,
    stage = STAGE,
    method = SAMPLING_METHOD,
    breeding_habitat = BREEDING_HABITAT
  ) %>% 
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  )


saveRDS(
  pa_data,
  "output/tabular/stephensi_pa_data.RDS"
)


# convert to time series of 

# col for ever detected
# col for earlierst detection or for last year if not detected

min_year <- min(pa_data$year, na.rm = TRUE) # 1984
max_year <- max(pa_data$year, na.rm = TRUE) # 2022

first_detection <- pa_data %>%
  select(x, y, year, presence) %>%
  mutate(
    year = if_else(is.na(year), min_year, year) # number of presen ces in India, Pakistan, Iran have no dates, assign as first year
  ) %>%
  group_by(x, y)  %>%
  summarise(
    presence = sum(presence),
    year = min(year),
    .groups = "drop"
  ) %>%
  mutate(
    year_first_detected = ifelse(presence == 0, 2022, year),
    ever_detected = ifelse(presence > 0, 1, 0)
  ) %>%
  select(
    x,
    y,
    year_first_detected,
    ever_detected 
  ) %>% 
  arrange(year_first_detected, ever_detected)

saveRDS(
  first_detection,
  "output/tabular/first_detection.RDS"
)

write_csv(
  first_detection,
  "output/tabular/first_detection.csv"
)
