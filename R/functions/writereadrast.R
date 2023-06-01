writereadrast <- function(x, filename){
  
  terra::writeRaster(
    x = x,
    filename = filename,
    overwrite = TRUE
  )
  
  terra::rast(x = filename)
  
}