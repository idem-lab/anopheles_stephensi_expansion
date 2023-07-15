library(readxl)
library(dplyr)

sapply(
  list.files("R/functions/", full.names = TRUE),
  source
)

# data from https://apps.who.int/malaria/maps/threats/
points_raw <- read_xlsx(
  path = "data/tabular/MTM_INVASIVE_VECTOR_SPECIES_20230607.xlsx",
  sheet = "Data"
)

aa <- points_raw %>%
  mutate(
    pa = case_when(
      is.na(MOSQUITO_NUMBER) ~ "Not 0",
      MOSQUITO_NUMBER == "0" ~ "0",
      TRUE ~ "Not 0"
    )
  )

table(aa$COUNTRY_NAME, aa$pa)

table(aa$SAMPLING_METHOD, aa$pa)

table(aa$INVASIVE_STATUS, aa$pa)

table(aa$COUNTRY_NAME, aa$INVASIVE_STATUS)

# next step dig into the 0s - should be plenty in Lao and Ethiopia, but not seeing any.
# Meanwhile can't se any of these 
