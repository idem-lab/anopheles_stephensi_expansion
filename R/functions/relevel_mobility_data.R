#' Takes a mobility data matrix, its associated location polygons in sf format,
#' and a desired new polygon set, relevels the mobility data to the level of the
#' new polygon set. The mobility data will be recalculated with weighting either
#' based on the area intersection with the original polygons, or based on
#' population proportion. That is, the original amount of mobility will be
#' distributed to new polygon i,j pairs, according to how much area or
#' population the new polygons share with original polygons.
#'
#' colnames are fixed for now, SORRY! will change
#'
#' @param original_geometry #expects an sf object
#' @param new_geometry  #expects an sf object
#' @param mobility_flow_matrix
#' @param weighting_by
#' @param pop_raster
#'
#' @return
#' @export
#'
#' @examples
relevel_mobility_data <- function(original_geometry, 
                                  #original_geometry_unique_id,
                                  new_geometry,
                                  #new_geometry_unique_id,
                                  mobility_flow_matrix,
                                  weighting_by = c("pop","area"),
                                  pop_raster = NULL) {
  match.arg(weighting_by)
  
  #first build intersection of the geometries, which is independent of weighting
  #might be a quirk with h3 to require this, not sure why
  #sf_use_s2(FALSE)
  
  #check if the column names match
  if (!("old_id" %in% names(original_geometry))) stop("column old_id not in original_geometry")
  if (!("new_id" %in% names(new_geometry))) stop("column old_id not in original_geometry")

  
  #intersect the two geometries
  intersect_geometry <- sf::st_intersection(new_geometry,original_geometry)
  
  if (weighting_by == "pop") {
    message("performing re-leleving with mobility data weighted by source and destination population")
    
    if (is.null(pop_raster)) stop("weighting by population but pop raster is null!")
    
    #use a pop raster to sample population in the insections
    intersect_zone_pop <- terra::zonal(pop_raster,
                             vect(intersect_geometry),
                             fun = sum,na.rm=TRUE)
    
    #add to intersect geometries 
    intersect_geometry$pop <- intersect_zone_pop[, 1]
    
    #get pops of original geometries
    intersect_pop_total <- intersect_geometry %>% #note we need use_pipe() when building this package!!!
      sf::st_drop_geometry() %>% 
      dplyr::group_by(old_id) %>% 
      dplyr::summarise(total.pop = sum(pop))
    
    intersect_geometry <- intersect_geometry %>% 
      sf::st_drop_geometry() %>% 
      dplyr::left_join(intersect_pop_total) %>% 
      #apply weight
      mutate(weight = pop/total.pop)

    
  } else if (weighting_by == "area") {
    message("performing re-leleving with mobility data weighted by source and destination area")
    
    #get area of interset
    intersect_geometry$area <- sf::st_area(intersect_geometry)
    
    #get areas of original geometries
    intersect_area_total <- intersect_geometry %>% 
      sf::st_drop_geometry() %>% 
      dplyr::group_by(old_id) %>% 
      dplyr::summarise(total.area = sum(area))
    
    intersect_geometry <- intersect_geometry %>% 
      sf::st_drop_geometry() %>% 
      dplyr::left_join(intersect_area_total) %>% 
      #apply weight
      mutate(weight = area/total.area)
    

  } else {
    stop("unknown weighting method!")
  }
  
  
  intersect_mat <- reshape2::dcast(intersect_geometry,
                                   old_id ~ new_id,
                                    value.var = "weight", 
                                    fun.aggregate = mean,
                                    fill=0)
  
  rownames(intersect_mat) <- intersect_mat[["old_id"]]
  intersect_mat <- as.matrix(intersect_mat[,2:ncol(intersect_mat)])
  
  #tally # of geos
  n_original_geometries <- nrow(intersect_mat)
  n_new_geometries <- ncol(intersect_mat)
  
  new_flow_mat <- matrix(0,
                         nrow = n_new_geometries, 
                         ncol = n_new_geometries)
  
  #iterate across original geometries to apply weighting
  for (i in 1:n_original_geometries) {
    for (j in 1:n_original_geometries) {
      
      this_flow_mat <- sapply(mobility_flow_matrix[i,j] * intersect_mat[i,],
                              FUN = function(x) x * intersect_mat[j,])
      
      new_flow_mat <- new_flow_mat + this_flow_mat
    }
  }
  
  
  return(new_flow_mat)
}
  