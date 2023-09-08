#devtools::install_github('SEEG-Oxford/movement')

library(tidyverse)
library(readr)
library(sf)
library(terra)
library(maptools)
library(movement)
library(geodata)
library(h3)
library(reshape2)

hexes <- readRDS("output/hexes.RDS")

#get wopr for all countries

global.pop <- rast("output/rasters/covariates/pop.grd")[["pop_2020"]]

global.pop <- terra::aggregate(global.pop,
                               fact = 10,
                               fun="sum",
                               na.rm = TRUE)

global.pop <- terra::project(global.pop,"ESRI:102022")

#make a spat vect version for calculating sum of pop
region_hex_spat_vec <- vect(hexes)


region_hex_pop <- terra::zonal(
  global.pop,
  region_hex_spat_vec,
  fun = sum,
  na.rm = TRUE)

#subsets to those with finite population count

# populated_hexes <- region_hex[!is.nan(region_hex_pop[,1]),]
# plot(populated_hexes)
# region_hex_pop <- region_hex_pop[!is.nan(region_hex_pop[,1]),]

#build a location data frame
location_data <- data.frame(location = hexes$h3_index,
                            population = region_hex_pop[, 1],
                            x = (st_centroid(hexes) %>% st_coordinates())[,1],
                            y = (st_centroid(hexes) %>% st_coordinates())[,2])

location_data  <- as.location_dataframe(location_data)

# simulate movements: note that the theta parameters (proportion of population
# moving) do not matter here since we are rescaling the movement matrix to
# have a maximum 1 anyway
predicted_flux_radiation  <- predict(originalRadiation(theta = 1),
                                     location_data,
                                     symmetric = FALSE)
radiation_matrix <- predicted_flux_radiation$movement_matrix

predicted_flux_gravity  <- predict(gravity(theta = 1,
                                           alpha = 1,
                                           beta = 1,
                                           gamma = 1),
                                   location_data,
                                   symmetric = FALSE)
gravity_matrix <- predicted_flux_gravity$movement_matrix

# do a straightforward distance matrix, for exponential distance decay too
distance_matrix <- fields::rdist.earth(location_data[, c("x", "y")])
rownames(distance_matrix) <- colnames(distance_matrix) <- rownames(predicted_flux_gravity$movement_matrix)

# scale this only against the maximum
radiation_matrix <- radiation_matrix / max(radiation_matrix)
gravity_matrix <- gravity_matrix / max(gravity_matrix)
distance_matrix <- distance_matrix / max(distance_matrix)

saveRDS(radiation_matrix, "output/tabular/radiation_matrix.RDS")
saveRDS(gravity_matrix, "output/tabular/gravity_matrix.RDS")
saveRDS(distance_matrix, "output/tabular/distance_matrix.RDS")
