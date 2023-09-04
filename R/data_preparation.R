library(dplyr)
library(readxl)
library(lubridate)


# data from https://apps.who.int/malaria/maps/threats
mtm_raw <- read_xlsx(
  path = "data/MTM_INVASIVE_VECTOR_SPECIES_20230303.xlsx",
  sheet = "Data"
)


pa_data <- mtm_raw %>%
  filter(VECTOR_SPECIES_COMPLEX == "An. stephensi" | VECTOR_SPECIES == "An. stephensi") %>% # the papers cited here have no mention of stephensi, unclear why in this data set
  mutate(
    #date = paste(YEAR_END, MONTH_END) %>% ym
    date = case_when(
      !is.na(MONTH_START) ~ paste(YEAR_START, MONTH_START) %>% ym,
      !is.na(YEAR_START) ~  paste(YEAR_START, "JAN") %>% ym,
      TRUE ~ NA_Date_
    ),
    presence = if_else(INVASIVE_STATUS == "not found", 0, 1) # either listed as "invasive", or "native" if present, or "not found" if absent
  ) %>%
  select(
    date,
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
