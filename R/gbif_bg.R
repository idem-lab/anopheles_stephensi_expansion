library(dplyr)
library(rgbif)
library(readr)
library(terra)
library(countrycode)

source("R/functions/maskpointsdf.R")

# based on bounding box code, which borrows from the nichemapR code,

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
    "KHM",
    "LAO",
    "MDV",
    "MMR",
    "NPL",
    "LKA",
    "THA",
    "TLS",
    "VNM"
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
#exclude the following countries
countries_exclude <- function() {
  c(
    "PRK",
    "IDN",
    "TLS",
    "EGY",
    "LBY",
    "MAR",
    "TUN"
  )
}


bg_countries <- function() {
  setdiff(countries(), countries_exclude())
}

# switch to 2-character codes for gbif
bg_countries_2 <- countrycode(
  bg_countries(),
  origin = "iso3c",
  destination = "iso2c"
)

tibble(
  country = countrycode(
    bg_countries(),
    origin = "iso3c",
    destination = "country.name"
  ),
  iso2 = bg_countries_2,
  iso3 = bg_countries()
) %>% print(n = 100)


# set up github credentials
# https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
# usethis::edit_r_environ()
# specify and save the below into the .Renviron file
# GBIF_USER="jwaller"
# GBIF_PWD="safe_fake_password_123"
# GBIF_EMAIL="jwaller@gbif.org"


# download animalia for focal region
# only needs to be run once, after this can just download per
# code below with occ_download_get

gbif_data  <- occ_download(
  pred('taxonKey', 1),
  pred_in('basisOfRecord',
          c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION")),
  pred_in('country', bg_countries_2),
  pred('hasGeospatialIssue', "FALSE"),
  pred('occurrenceStatus', "PRESENT"),
  pred("hasCoordinate", TRUE),
  pred_lt("coordinateUncertaintyInMeters",1000),
  pred_gte('year', 2010),
  format = "SIMPLE_CSV"
)

gbif_data
# 
# occ_download_wait('0005642-230828120925497', curlopts=list(http_version=2))
# re curlopts seems to be necessary possibly only on mac:
# https://github.com/ropensci/rgbif/issues/579


gbif_citation("0005642-230828120925497")

bg_ani_gbif_raw <- occ_download_get('0005642-230828120925497') %>%
  occ_download_import() %>% 
  write_csv(
    sprintf(
      "data/tabular/gbif_region_animalia_%s.csv",
      lubridate::today() %>%
        format("%Y%m%d")
    )
  )


bg_ani_df_raw <- bg_ani_gbif_raw %>%
  filter(species != "Anopheles stephensi") %>%
  dplyr::select(
    lon = decimalLongitude,
    lat = decimalLatitude
  ) %>%
    distinct


covmask <- rast(x = "output/rasters/covariates/covmask.grd")

bg_ani <- maskpointsdf(
  df = bg_ani_df_raw,
  msk = covmask
)

plot(covmask)
points(vect(bg_ani))

saveRDS(
  bg_ani,
  sprintf(
    "output/tabular/bg_animalia_%s.RDS",
    lubridate::today() %>%
      format("%Y%m%d")
  )
)

#### anopheles aedes culex only


gbif_moz  <- occ_download(
  pred_or(
    pred("taxonKey", 7924646), # Aedes
    pred("taxonKey", 1650098), # Anopheles
    pred("taxonKey", 1497010) # Culex
  ),
  pred_in('country', bg_countries_2),
  pred("hasCoordinate", TRUE),
  format = "SIMPLE_CSV"
)

gbif_moz

occ_download_wait('0005643-230828120925497', curlopts=list(http_version=2))

gbif_citation("0005643-230828120925497")

bg_moz_gbif_raw <- occ_download_get('0005643-230828120925497') %>%
  occ_download_import() %>% 
  write_csv(
    sprintf(
      "data/tabular/gbif_region_moz_%s.csv",
      lubridate::today() %>%
        format("%Y%m%d")
    )
  )


# subset to only results with low coordinate uncertainty (from all years)
bg_moz_best_raw <- bg_moz_gbif_raw %>%
  filter(coordinateUncertaintyInMeters < 1000) %>%
  filter(species != "Anopheles stephensi") %>%
  dplyr::select(
    lat = decimalLatitude,
    lon = decimalLongitude
  ) %>%
  distinct

# subset to data post 2000 with either low coordinate uncertainty or missing
# coordinate uncertainty (logic here is that most modern points will come from GPS
# so should be tolerably accurate unless specified otherwise)
bg_moz_many_raw <- bg_moz_gbif_raw %>%
  filter(year >= 2000) %>%
  filter(
    is.na(coordinateUncertaintyInMeters)
    | coordinateUncertaintyInMeters < 1000
  ) %>%
  filter(species != "Anopheles stephensi") %>%
  dplyr::select(
    lat = decimalLatitude,
    lon = decimalLongitude
  ) %>%
  distinct


# mask to land in area of interest
bg_moz_best <- maskpointsdf(
  df = bg_moz_best_raw,
  msk = covmask
)

bg_moz_many <- maskpointsdf(
  df = bg_moz_many_raw,
  msk = covmask
)

plot(covmask)
points(vect(bg_ani))
points(vect(bg_moz_many), col = "grey60")
points(vect(bg_moz_best), col = "orange")

# save
saveRDS(
  bg_moz_best,
  sprintf(
    "output/tabular/bg_moz_best_%s.RDS",
    lubridate::today() %>%
      format("%Y%m%d")
  )
)

saveRDS(
  bg_moz_many,
  sprintf(
    "output/tabular/bg_moz_many_%s.RDS",
    lubridate::today() %>%
      format("%Y%m%d")
  )
)


