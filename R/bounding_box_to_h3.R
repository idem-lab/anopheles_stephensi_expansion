source("R/bounding_box.R")

lapply(
  list.files("R/functions"),
  source
)

region_hex <- get_h3_from_sf(
  sf_object = st_as_sf(region_shape_buffer))

plot(region_hex)

