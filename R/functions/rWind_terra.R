#rWind functions converted from raster to terra

wind2terra <- function(x){
  if (inherits(x, "rWind_series")) {
    X <- lapply(x, wind2terra_int)
  }
  else {
    return(wind2terra_int(x))
  }
  X
}

wind2terra_int <- function(x){
  ras_dir <- rast(
    x[, c("lon", "lat", "dir")] %>%
      as.data.frame,
    type = "xyz",
    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  )
  ras_speed <- rast(
    x[, c("lon", "lat", "speed")] %>%
      as.data.frame,
    type = "xyz",
    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  )
  
  tmp <- c(ras_dir, ras_speed)
  names(tmp) <- c("direction", "speed")
  return(tmp)
}



flow.dispersion <- function (x, fun = cost.FMGS, output = "transitionLayer", ...) 
{
  flow.dispersion_int(x, fun = fun, output = output, 
                               ...)
  
}


flow.dispersion_int_terra <- function(stack, fun = cost.FMGS, output = "transitionLayer",
                                ...) {
  output <- match.arg(output, c("raw", "transitionLayer"))
  
  DL <- as.matrix(stack$direction)
  SL <- as.matrix(stack$speed)
  M <- matrix(as.integer(1:ncell(stack$direction)),
              nrow = nrow(stack$direction), byrow = TRUE
  )
  nr <- nrow(M)
  nc <- ncol(M)
  
  ###################################################################
  
  directions <- c(315, 0, 45, 270, 90, 225, 180, 135)
  
  ###################################################################
  
  # Go Nortwest
  
  north.west.from <- as.vector(M[-1, -1])
  north.west.to <- as.vector(M[-nr, -nc])
  north.west.cost <- fun(DL[-1, -1], SL[-1, -1], directions[1], ...)
  
  ###################################################################
  
  # Go North
  
  north.from <- as.vector(M[-1, ])
  north.to <- as.vector(M[-nr, ])
  north.cost <- as.vector(fun(DL[-1, ], SL[-1, ], directions[2], ...))
  
  ###################################################################
  
  # Go Norteast
  
  north.east.from <- as.vector(M[-1, -nc])
  north.east.to <- as.vector(M[-nr, -1])
  north.east.cost <- as.vector(fun(DL[-1, -nc], SL[-1, -nc], directions[3], ...))
  
  ###################################################################
  
  # Go West
  
  west.from <- as.vector(M[, -1])
  west.to <- as.vector(M[, -nc])
  west.cost <- as.vector(fun(DL[, -1], SL[, -1], directions[4], ...))
  
  ###################################################################
  
  # Go East
  
  east.from <- as.vector(M[, -nc])
  east.to <- as.vector(M[, -1])
  east.cost <- as.vector(fun(DL[, -nc], SL[, -nc], directions[5], ...))
  
  ###################################################################
  
  # Go Southwest
  
  south.west.from <- as.vector(M[-nr, -1])
  south.west.to <- as.vector(M[-1, -nc])
  south.west.cost <- as.vector(fun(DL[-nr, -1], SL[-nr, -1], directions[6], ...))
  
  ###################################################################
  
  # Go South
  
  south.from <- as.vector(M[-nr, ])
  south.to <- as.vector(M[-1, ])
  south.cost <- as.vector(fun(DL[-nr, ], SL[-nr, ], directions[7], ...))
  
  ###################################################################
  
  # Go Southeast
  
  south.east.from <- as.vector(M[-nr, -nc])
  south.east.to <- as.vector(M[-1, -1])
  south.east.cost <- as.vector(fun(DL[-nr, -nc], SL[-nr, -nc], directions[8], ...))
  
  ###################################################################
  
  ii <- c(north.west.from, north.from, north.east.from, west.from, east.from, south.west.from, south.from, south.east.from)
  jj <- c(north.west.to, north.to, north.east.to, west.to, east.to, south.west.to, south.to, south.east.to)
  xx <- c(north.west.cost, north.cost, north.east.cost, west.cost, east.cost, south.west.cost, south.cost, south.east.cost)
  
  tl <- sparseMatrix(i = ii, j = jj, x = xx)
  if (output == "raw") {
    return(tl)
  }
  if (output == "transitionLayer") {
    tmp <- raster::transition(raster::raster(stack$direction), transitionFunction = function(x) 0, directions = 8)
    transitionMatrix(tmp) <- sparseMatrix(i = ii, j = jj, x = 1 / xx)
    return(tmp)
  }
  return(NULL)
}


flow.dispersion_int_terra <- function(stack, fun = cost.FMGS, output = "transitionLayer",
                                      ...) {
  output <- match.arg(output, c("raw", "transitionLayer"))
  
  DL <- as.matrix(stack$direction)
  SL <- as.matrix(stack$speed)
  M <- matrix(as.integer(1:ncell(stack$direction)),
              nrow = nrow(stack$direction), byrow = TRUE
  )
  nr <- nrow(M)
  nc <- ncol(M)
  
  ###################################################################
  
  directions <- c(315, 0, 45, 270, 90, 225, 180, 135)
  
  ###################################################################
  
  # Go Nortwest
  
  north.west.from <- as.vector(M[-1, -1])
  north.west.to <- as.vector(M[-nr, -nc])
  north.west.cost <- fun(DL[-1, -1], SL[-1, -1], directions[1])
  
  ###################################################################
  
  # Go North
  
  north.from <- as.vector(M[-1, ])
  north.to <- as.vector(M[-nr, ])
  north.cost <- as.vector(fun(DL[-1, ], SL[-1, ], directions[2]))
  
  ###################################################################
  
  # Go Norteast
  
  north.east.from <- as.vector(M[-1, -nc])
  north.east.to <- as.vector(M[-nr, -1])
  north.east.cost <- as.vector(fun(DL[-1, -nc], SL[-1, -nc], directions[3]))
  
  ###################################################################
  
  # Go West
  
  west.from <- as.vector(M[, -1])
  west.to <- as.vector(M[, -nc])
  west.cost <- as.vector(fun(DL[, -1], SL[, -1], directions[4]))
  
  ###################################################################
  
  # Go East
  
  east.from <- as.vector(M[, -nc])
  east.to <- as.vector(M[, -1])
  east.cost <- as.vector(fun(DL[, -nc], SL[, -nc], directions[5]))
  
  ###################################################################
  
  # Go Southwest
  
  south.west.from <- as.vector(M[-nr, -1])
  south.west.to <- as.vector(M[-1, -nc])
  south.west.cost <- as.vector(fun(DL[-nr, -1], SL[-nr, -1], directions[6]))
  
  ###################################################################
  
  # Go South
  
  south.from <- as.vector(M[-nr, ])
  south.to <- as.vector(M[-1, ])
  south.cost <- as.vector(fun(DL[-nr, ], SL[-nr, ], directions[7]))
  
  ###################################################################
  
  # Go Southeast
  
  south.east.from <- as.vector(M[-nr, -nc])
  south.east.to <- as.vector(M[-1, -1])
  south.east.cost <- as.vector(fun(DL[-nr, -nc], SL[-nr, -nc], directions[8]))
  
  ###################################################################
  
  ii <- c(north.west.from, north.from, north.east.from, west.from, east.from, south.west.from, south.from, south.east.from)
  jj <- c(north.west.to, north.to, north.east.to, west.to, east.to, south.west.to, south.to, south.east.to)
  xx <- c(north.west.cost, north.cost, north.east.cost, west.cost, east.cost, south.west.cost, south.cost, south.east.cost)
  
  tl <- sparseMatrix(i = ii, j = jj, x = xx)
  if (output == "raw") {
    return(tl)
  }
  if (output == "transitionLayer") {
    tmp <- raster::transition(raster::raster(stack$direction), transitionFunction = function(x) 0, directions = 8)
    transitionMatrix(tmp) <- sparseMatrix(i = ii, j = jj, x = 1 / xx)
    return(tmp)
  }
  return(NULL)
}