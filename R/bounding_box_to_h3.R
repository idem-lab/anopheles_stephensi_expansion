#' takes a SF object, get corresponding h3 hexagons in pre-specified projection and resolution, as
#' SF polygons
#'
#' @param sf_object
#' @param h3_crs_to
#' @param h3_res
#'
#' @return
#' @export
#'
#' @examples
get_h3_from_sf <- function(
    sf_object,
    h3_crs_to = "ESRI:102022",
    h3_res = 1) {
  
  #get hexes ID filling up the desired polygon
  region_hex <- h3::polyfill(sf_object,
                             res = h3_res)
  
  #remove possible duplicates for hexes
  region_hex <- unique(region_hex)
  
  #get hexes in sf polygons
  region_hex <- h3::h3_to_geo_boundary_sf(region_hex)
  
  # fix hex poly projections
  region_hex <- sf::st_transform(region_hex,crs = h3_projection)
  
  return(region_hex)
  
}

source("R/bounding_box.R")
region_hex <- get_h3_from_sf(
  sf_object = st_as_sf(region_shape_buffer))

plot(region_hex)

