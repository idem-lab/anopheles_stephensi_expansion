library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(countrycode)


# africa countries with A stephensi older records
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


ve_raw <- read_csv("data/tabular/vector_extraction_data.csv") %>%
  dplyr::select(
    month_start = sample_period_month_start,
    year_start = sample_period_year_start,
    month_end = sample_period_month_end,
    year_end = sample_period_year_end,
   count = collection_count,
   method = collection_method_name,
   sp = anopheline_species,
   country = vector_site_country,
   lon = vector_site_coordinates_longitude,
   lat = vector_site_coordinates_latitude
  )


seasonality <- ve_raw %>%
  drop_na %>%
  filter(
    year_start %in% 1900:2015,
    year_end %in% 1900:2015,
    month_start %in% 1:12,
    month_end %in% 1:12,
  ) %>%
  mutate(
    start = sprintf(
      "%s%02d%02d",
      year_start,
      month_start,
      1
    ) %>%
      ymd,
    end = sprintf(
      "%s%02d%02d",
      year_end,
      month_end,
      1
    ) %>%
      ymd +
      months(1) -
      days(1),
    survey_length = difftime(end, start, units = "days")
  ) %>% 
  filter(survey_length <= 93) %>%
  mutate(
    date = start + survey_length/2,
    month = month(date)
  ) %>%
  filter(
    country %in% countrycode(
      sourcevar = africa_countries(),
      origin = "iso3c",
      destination = "country.name"
    )
  )



ggplot(seasonality) +
  geom_point(
    aes(
      x = month,
      y = count,
      colour = country
    )
  )

ggplot(
  seasonality
) +
  geom_point(
    aes(
      x = month,
      y = count,
      colour = country
    )
  ) +
  geom_smooth(
    aes(
      x = month,
      y = count,
      colour = country
    ),
    method = "loess"
  ) +
  facet_wrap(
    ~ country,
    scales = "free"
  )

