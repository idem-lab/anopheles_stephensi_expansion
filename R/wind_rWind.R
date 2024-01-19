library(rWind)
library(gdistance)
library(dplyr)
library(terra)
library(tibble)
library(tidyr)

#data("wind.data")
#wind.data <- as_tibble(wind.data)


source("R/functions/rWind_terra.R")

# download wind raster for date / time
w <- wind.dl(2023, 2, 12, 12, -7, -4, 34.5, 37.5)

# convert to direction and speed from u and v components
wind_layer <- wind2raster(w)

# convert to friction layer
conductance <- flow.dispersion(wind_layer)

# create long list of all cells vs all cells and get cost from
# movement from every cell to every other
locs <- w %>%
  select(lon, lat) %>%
  as_tibble

locs1 <- locs %>%
  rename(lon1 = lon, lat1 = lat)
locs2 <- locs %>%
  rename(lon2 = lon, lat2 = lat)

library(purrr)

z <- function(x,y){data.frame(x,y)}

expand_grid(locs1, locs2) %>%
  mutate(
    l1 = mapply(c, lon1, lat1, SIMPLIFY = FALSE),
    l2 = mapply(c, lon2, lat2, SIMPLIFY = FALSE),
    cost = mapply(
      FUN = costDistance,
      fromCoords = l1,
      toCoords = l2,
      MoreArgs = list(x = conductance),
      SIMPLIFY = FALSE
    )
  ) %>% 
  unnest(cost)
  

