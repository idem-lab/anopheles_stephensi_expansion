



relevel_mobility_data <- function(original_geometry, 
                                  new_geometry,
                                  weighting_by = c("pop","area"),
                                  pop_raster = NULL) {
  match.arg(weighting_by)
  
  if (weighting_by == "pop") {
    message("performing re-leleving with mobility data weighted by source and destination population")
    
    if (is.null(pop_raster)) stop("weighting by population but pop raster is null!")
    #use a pop raster to sample population in the original catchments
    origin_zone_pop <- zonal(pop_raster,
                             original_geometry,
                             fun = sum,na.rm=TRUE)
    
    #add to original zones
    original_geometry$population <- origin_zone_pop
    
  } else if (weighting_by == "area") {
    message("performing re-leleving with mobility data weighted by source and destination area")
  } else {
    stop("unknown weighting method!")
  }
  
}
  


library(tidyverse)
library(readr); library(sf); library(terra); library(maptools); library(movement); library(geodata); library(h3);library(reshape2)

#get wopr for all countries

global.pop <- rast("data/ppp_2020_1km_Aggregated.tif")
plot(global.pop)

#get flow data
international.migration.flow.data <- read_csv("data/SexDisaggregated_Migration/MigrationEstimates/MigrEst_international_v7.csv")
subnational.migration.flow.data <- read_csv("data/SexDisaggregated_Migration/MigrationEstimates/MigrEst_internal_v4.csv")

#target ext
target.ext <- ext(-32.9770992366411, 100.30534351145, -36.7374249758784, 43.4248387475794)

#get centroids
international.migration.flow.geometry <- vect("data/SexDisaggregated_Migration/SpatialData/")
international.migration.flow.geometry <- international.migration.flow.geometry %>% 
  subset(international.migration.flow.geometry$CONT %in% c("AFR","ASI"))


international.migration.flow.geometry <- crop(international.migration.flow.geometry,target.ext)

plot(international.migration.flow.geometry)

#get country lists
countries <- unique(international.migration.flow.geometry$ISO)

countries.with.ssd <- countries
countries.with.ssd[countries.with.ssd=="SUD"] <- "SDN"
countries.with.ssd <- c(countries.with.ssd,"SSD","ESH")
#get admin1 polygon


admin1.outline <- gadm(countries.with.ssd, level=1, path = "data/gadm/", version="latest", resolution=1)



#   vect("data/ne_10m_admin_1_states_provinces/")
# admin1.outline <- admin1.outline[,c("adm0_a3","adm1_code")]
# 
# admin1.outline <- admin1.outline %>% subset(admin1.outline$adm0_a3 %in% countries)

admin1.outline <- crop(admin1.outline,target.ext)

#plot(admin1.outline)


admin1.pop <- zonal(global.pop,
                    admin1.outline,
                    fun = sum,na.rm=TRUE)

admin1.outline$population <- admin1.pop
admin1.outline$weight <- admin1.outline$population / sum(admin1.outline$population, na.rm = TRUE)

#plot(admin1.outline,"weight")


# 
# flowlines <- vect("data/SexDisaggregated_Migration/SpatialData/Flowlines/international/")
# flowlines <- crop(flowlines,target.ext)
# plot(flowlines)

test <- st_as_sf(international.migration.flow.geometry)

test.hex <- geo_to_h3(test, res = 0)
test.hex <- unique(test.hex)
test.hex <- h3_to_geo_boundary_sf(test.hex)
plot(test.hex)



admin0.outline <- gadm(countries.with.ssd, level=0, 
                       path = "data/gadm/", 
                       version="latest", 
                       resolution=1)

admin0.outline <- crop(admin0.outline,target.ext)

#plot(admin0.outline)

test <- st_as_sf(admin0.outline)

test <- test[unique(test$COUNTRY),]
#plot(test)
sf_use_s2(FALSE)
# test <- st_make_valid(test)
# test.hex <- st_make_valid(test.hex)

test_intersect <- st_intersection(test.hex,test)

plot(test_intersect[1:10,])


test_intersect %>% filter(COUNTRY == "Angola") %>% summarise(n = length(unique(h3_index)))

test_intersect$area <- st_area(test_intersect)

test_intersect.total <- test_intersect %>% as_tibble() %>% 
  group_by(COUNTRY) %>% 
  summarise(total.area = sum(area))


test_intersect <- test_intersect %>% 
  left_join(test_intersect.total)

test_intersect$weight <- test_intersect$area/test_intersect$total.area
test_intersect$weight <- as.numeric(test_intersect$weight)


test_intersect_tibble <- test_intersect %>% st_drop_geometry()

test_intersect_tibble.1 <- dcast(test_intersect_tibble,
                                 COUNTRY ~ h3_index,
                                 value.var = "weight", 
                                 fun.aggregate = mean,
                                 fill=0)

rownames(test_intersect_tibble.1) <- test_intersect_tibble.1$COUNTRY
intersect_mat <- as.matrix(test_intersect_tibble.1[,2:26])

fake_flow_mat <- matrix(1,82,82)

fake_flow_mat * intersect_mat %>% View
#test_intersect %>% as_tibble() %>% group_by(COUNTRY) %>% summarise(sum(weight))

test.hex <- polyfill(country.outline, res = 1)
test.hex <- h3_to_geo_boundary_sf(test.hex)
plot(test.hex)
