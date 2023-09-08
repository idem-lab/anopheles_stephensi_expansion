#devtools::install_github('SEEG-Oxford/movement')

library(tidyverse)
library(readr); library(sf); library(terra); library(maptools); library(movement); library(geodata); library(h3);library(reshape2)

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
# simulate movements 
#note that the theta parameter (proportion of population moving) does not matter here since we are rescaling every column of the movement matrix to sum to 1 anyway
predicted_flux  <- predict(radiationWithSelection(), location_data, symmetric = FALSE)

# scale this only against the maximum
movement_matrix <- predicted_flux$movement_matrix / max(predicted_flux$movement_matrix)

# zero-out the diagonals

#check if sums are sensible 
colSums(movement_matrix)
# # fit a new model to these data
# movement_model <- movement(movement_matrix ~ location_data, radiationWithSelection())
# # print movement_model
# print(movement_model)
# # predict the population movements
# predicted_movements  <- predict(movement_model, location_data)
# # display the predicted movements
# plot(predicted_movements)

saveRDS(movement_matrix, "output/tabular/radiation_matrix.RDS")

