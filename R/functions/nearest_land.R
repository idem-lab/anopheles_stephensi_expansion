#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param coords
#' @param raster
#' @param projected_crs
#' @param max_distance_m
#' @return
#' @author Nick Golding
#' @export
nearest_land <- function(coords,
                         raster,
                         projected_crs = "ESRI:102022",
                         max_distance_m = 10000) {
  
  old_crs <- crs(raster)
  
  # find the index to cells that need to be relocated
  coords_new <- coords %>%
    st_as_sf(
      crs = old_crs,
      coords = c("x", "y")
    ) %>%
    mutate(
      point = row_number(),
      valid = !is.na(terra::extract(raster, .)[, 2]),
      .before = everything()
    )
  
  if (all(coords_new$valid)) {
    return(coords_new)
  }
  
  # project the raster
  raster_proj <- raster %>%
    terra::project(projected_crs) %>%
    `names<-`("raster")
  
  # project and buffer the coordinates that are not in valid cells
  coords_proj <- coords_new %>%
    filter(
      !valid
    ) %>%
    st_transform(
      projected_crs
    )
  
  buffered <- coords_proj %>%
    st_buffer(max_distance_m)
  
  # find all the points within that buffer, and convert to coordinates (in old CRS)
  cells <- terra::extract(raster_proj,
                          buffered,
                          cells = TRUE) %>%
    bind_cols(
      terra::xyFromCell(raster_proj, .$cell)
    ) %>%
    st_as_sf(
      crs = projected_crs,
      coords = c("x", "y")
    ) %>%
    filter(
      !is.na(raster)
    ) %>%
    bind_cols(
      st_geometry(.) %>%
        st_transform(old_crs) %>%
        st_coordinates()
    ) %>%
    # join on the original point number (ID is the index to the subset)
    left_join(
      tibble(
        point = coords_proj$point,
        id = seq_len(nrow(coords_proj))
      ),
      by = c(ID = "id")
    ) %>%
    dplyr::select(
      point,
      x_proposed = X,
      y_proposed = Y
    ) %>%
    st_drop_geometry()
  
  # return the nearest point on land (or NA if none could be found)
  coords_new %>%
    mutate(
      x = st_coordinates(.)[, 1],
      y = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry() %>%
    left_join(
      cells,
      by = "point",
      multiple = "all"
    ) %>%
    rowwise() %>%
    mutate(
      distance = fields::rdist.earth(
        cbind(x, y),
        cbind(x_proposed, y_proposed)
      )[, 1]
    ) %>%
    group_by(point, valid, x, y) %>%
    mutate(
      keep = case_when(
        # don't need a new point
        valid ~ TRUE,
        # need a new point and found one
        !valid & distance == min(distance) ~ TRUE,
        # need but didn't find a new point (take first NA)
        !valid & all(is.na(distance)) & row_number() == 1 ~ TRUE,
        .default = FALSE
      )
    ) %>%
    filter(keep) %>%
    mutate(
      x = case_when(
        valid ~ x,
        !valid & keep ~ x_proposed,
        .default = NA
      ),
      y = case_when(
        valid ~ y,
        !valid & keep ~ y_proposed,
        .default = NA
      )
    ) %>%
    ungroup() %>%
    dplyr::select(
      x,
      y
    )
  
}
