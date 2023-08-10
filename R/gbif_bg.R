library(dplyr)
library(rgbif)
library(readr)

# based on bounding box code, which borrows from the nichemapR code,

# africa countries with A stephensi older records
africa_countries <- function () {
  c(
    # "AGO",
    # "BDI",
    # "BEN",
    # "BFA",
    # "BWA",
    # "CAF",
    # "CIV",
    # "CMR",
    # "COD",
    # "COG",
    # "COM",
    # "CPV",
    "DJI",
    # "DZA",
    # "EGY",
    # "ERI",
    # "ESH",
    "ETH",
    # "GAB",
    # "GHA",
    # "GIN",
    # "GMB",
    # "GNB",
    # "GNQ",
    # "KEN",
    # "LBR",
    # "LBY",
    # "LSO",
    # "MAR",
    # "MDG",
    # "MLI",
    # "MOZ",
    # "MRT",
    # "MUS",
    # "MWI",
    # "NAM",
    # "NER",
    # "NGA",
    # "RWA",
    # "SDN",
    # "SEN",
    # "SLE",
    # "SOM",
    "SSD"#,
    # "STP",
    # "SWZ",
    # "TCD",
    # "TGO",
    # "TUN",
    # "TZA",
    # "UGA",
    # "ZAF",
    # "ZMB",
    # "ZWE"
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


bg_countries <- function() {
  setdiff(countries(), countries_exclude())
}

# switch to 2-character codes for gbif
bg_countries_2 <- countrycode::countrycode(
  bg_countries(),
  origin = "iso3c",
  destination = "iso2c"
)



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

# gbif_data  <- occ_download(
#   pred('taxonKey', 1),
#   pred_in('basisOfRecord', 
#           c("MACHINE_OBSERVATION", "HUMAN_OBSERVATION")),
#   pred_in('country', bg_countries_2),
#   pred('hasGeospatialIssue', "FALSE"),
#   pred('occurrenceStatus', "PRESENT"),
#   pred("hasCoordinate", TRUE),
#   pred_lt("coordinateUncertaintyInMeters",1000),
#   pred_gte('year', 2010),
#   format = "SIMPLE_CSV"
# )
# 
# gbif_data

gbif_citation("0013949-230530130749713")

gbif_bg_raw <- occ_download_get("0013949-230530130749713") %>%
  occ_download_import() %>% 
  write_csv(
    sprintf(
      "data/tabular/gbif_region_animalia_%s.csv",
      "20230607"
    )
  )


bg_points <- gbif_bg_raw %>%
  select(
    lat = decimalLatitude,
    lon = decimalLongitude
  )

saveRDS(
  bg_points,
  "output/bg_points.RDS"
)
