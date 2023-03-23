# given a raster 'pop_raster' where each cell value contains the population and
# a positive integer number of patches 'n_patches' required, return a list
# containing a SpatVector object 'patch_centroids' giving the patch centroid
# coorindates, IDs, and total populations, and a raster 'patch_raster' where
# values give the ID of the atch to which they belong. n_representative points
# is the number of random population-weighted points to use in the clustering.
patchify <- function(pop_raster, n_patches, n_representative_points = 2000) {
  
  require(terra)
  require(dplyr)
  
  # sample random points
  random_points <- terra::spatSample(pop_raster,
                                     n_representative_points,
                                     method = "weights",
                                     replace = TRUE,
                                     na.rm = TRUE,
                                     as.points = TRUE)
  
  # do k-means clustering on these
  kmn <- random_points %>%
    sf::st_as_sf() %>%
    sf::st_coordinates() %>%
    kmeans(centers = n_patches)
  
  # get the patch centroids and give them IDs and weights
  patch_centroids <- kmn$centers %>%
    vect(
      atts = data.frame(id = seq_len(n_patches))
    )
  
  # get a voronoi tesselation of this, on a raster
  patch_raster <- patch_centroids %>%
    voronoi() %>%
    rasterize(pop_raster,
              field = "id",
              fun = "min") %>%
    mask(pop_raster)
  
  # get population counts
  pop_totals <- terra::zonal(pop_raster, patch_raster, fun = sum)
  
  # this is dangerous - how do we join on the ID column?
  patch_centroids$population <- pop_totals[, 2]
  patch_centroids$weight <- patch_centroids$population / sum(patch_centroids$population)
  
  # return both things
  list(
    patch_centroids = patch_centroids,
    patch_raster = patch_raster
  )
}

# demo:

# get WorldPop population raster for Nigeria
# remotes::install_github('wpgp/wopr')
library(wopr)
catalogue <- getCatalogue()
nigeria_pop_ref <- subset(catalogue,
       country == 'NGA' &
         category == 'Population' &
         version == 'v1.2' & 
         filetype == 'gridded')

downloadData(nigeria_pop_ref, wopr_dir = "data")
unzip("data/NGA/population/v1.2/NGA_population_v1_2_gridded.zip")

library(terra)
nga_pop_100 <- rast("wopr/NGA/population/v1.2/NGA_population_v1_2_gridded/NGA_population_v1_2_gridded.tif")

# this is at 100m resolution, aggregate up to 5km resolution
nga_pop <- terra::aggregate(nga_pop_100, 50, fun = "sum", na.rm = TRUE)

patches <- patchify(nga_pop, 12)

# plot the results
par(mfrow = c(2, 1))
plot(nga_pop)
points(patches$patch_centroids,
       cex = 10 * patches$patch_centroids$weight)

plot(patches$patch_raster)
points(patches$patch_centroids,
       cex = 10 * patches$patch_centroids$weight)
