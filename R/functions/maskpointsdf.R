maskpointsdf <- function(df, msk){
  
  vctraw <- terra::vect(
    x = df,
    crs = "+proj=longlat +datum=WGS84"
  )
  
  pol <- terra::as.polygons(msk)
  
  vct <- terra::mask(vctraw, pol)
  
  terra::geom(vct) %>%
    as_tibble %>%
    select(
      lon = x,
      lat = y
    )
  
}