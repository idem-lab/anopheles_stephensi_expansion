library(remotes)
install_github("matthewkling/windscape")


### generic example code

library(windscape)
library(rWind)
library(raster)
library(tidyverse)
library(gdistance)
library(geosphere)

# download some wind time series data using the rWind package
# we'll grab 24 days distributed across a single year
times <- expand.grid(yyyy = 2023, mm = 1:12, dd = c(1, 15), tt = 12)
wd <- pmap(times, wind.dl, lon1 = -120, lon2 = -90, lat1 = 30, lat2 = 50)

# convert to raster stack
u <- wd %>%
  map(dplyr::select, lon, lat, ugrd10m) %>%
  map(rasterFromXYZ) %>%
  stack()
v <- wd %>%
  map(dplyr::select, lon, lat, vgrd10m) %>%
  map(rasterFromXYZ) %>%
  stack()
wind <- stack(u, v)

# convert wind time series into a wind conductance raster
conductance <- wind %>%
  windrose_rasters(outfile="windrose.tif", order="uuvv", p=1, ncores = 1) %>% 
  add_coords()


# fails 
# Error in .calcTest(x[1:5], fun, na.rm, forcefun, forceapply) : 
#   cannot use this function

# within windrose_rasters()

w <- wind
outfile <- "windrose.tif"
order <- "uuvv"
p <- 1
ncores <- 1


rosefun <- function(x) windrose(x, ...)
wr <- add_res(w)
wr <- add_lat(wr)
wr <- raster::calc(wr, fun = rosefun, forceapply = TRUE, 
                     filename = outfile)

# fails
# Error in .calcTest(x[1:5], fun, na.rm, forcefun, forceapply) : 
#   cannot use this function

# within windrose()
x <- wr

lat <- x[1]
res <- x[2]
x <- x[3:length(x)]
uv <- matrix(x, ncol = 2, byrow = F)
weight <- sqrt(uv[, 1]^2 + uv[, 2]^2)^p

zz <- windscape::direction(uv[, 2], -1 * uv[, 1])

dir <- spin90(windscape::direction(uv[, 2], -1 * uv[, 1]))

# fails
# Error in x[x < (-180)] <- x[x < (-180)] + 360 : 
#   NAs are not allowed in subscripted assignments

# alter spin90 to cope with NAs
spin90 <- function (x) 
{
  x <- x - 90
  
  ss <- x < (-180) & !is.na(x)
  
  x[ss] <- x[ss] + 360
  x
}

dir <- spin90(windscape::direction(uv[, 2], -1 * uv[, 1]))


dirgr0 <- dir < 0 & !is.na(dir)

dir[dirgr0] <- dir[dirgr0] + 360

dir0 <- dir == 0

dir[dir0] <- 360

nc <- cbind(x = c(0, res, res, res, 0, -res, -res, -res), 
            y = c(res, res, 0, -res, -res, -res, 0, res) + lat)
nc[, 2] <- pmin(nc[, 2], 90)
nc[, 2] <- pmax(nc[, 2], -90)
nb <- geosphere::bearingRhumb(c(0, lat), nc)
nb <- c(nb, 360)
l <- edge_loadings(dir, weight, nb)
nd <- geosphere::distGeo(c(0, lat), nc)
l <- l/nd
l[c(6:8, 1:5)]
}


x <- zz

spin90 <- function (x) 
{
  x <- x - 90
  
  ss <- x < (-180) & !is.na(x)
  
  x[ss] <- x[ss] + 360
  x
}


windrose <- function (x, p = 1) 
{
  lat <- x[1]
  res <- x[2]
  x <- x[3:length(x)]
  uv <- matrix(x, ncol = 2, byrow = F)
  weight <- sqrt(uv[, 1]^2 + uv[, 2]^2)^p
  dir <- spin90(windscape::direction(uv[, 2], -1 * uv[, 1]))
  dirgr0 <- dir < 0 & !is.na(dir)
  dir[dirgr0] <- dir[dirgr0] + 360
  dir0 <- dir == 0
  dir[dir0] <- 360
  nc <- cbind(x = c(0, res, res, res, 0, -res, -res, -res), 
              y = c(res, res, 0, -res, -res, -res, 0, res) + lat)
  nc[, 2] <- pmin(nc[, 2], 90)
  nc[, 2] <- pmax(nc[, 2], -90)
  nb <- geosphere::bearingRhumb(c(0, lat), nc)
  nb <- c(nb, 360)
  l <- edge_loadings(dir, weight, nb)
  nd <- geosphere::distGeo(c(0, lat), nc)
  l <- l/nd
  l[c(6:8, 1:5)]
}


function (w, outfile, order = "uvuv", ncores = 1, ...) 
{
  if (order == "uvuv") {
    even <- function(x) x%%2 == 0
    w <- w[[c(which(!even(1:nlayers(w))), which(even(1:nlayers(w))))]]
  }
  rosefun <- function(x) windrose(x, ...)
  if (ncores == 1) {
    wr <- add_res(w)
    wr <- add_lat(wr)
    wr <- raster::calc(wr, fun = rosefun, forceapply = TRUE, 
                       filename = outfile)
    names(wr) <- c("SW", "W", "NW", "N", "NE", "E", "SE", 
                   "S")
    return(wr)
  }
  else {
    nlayer <- nlayers(w)/2
    core <- rep(1:ncores, each = floor(nlayer/ncores))
    core <- c(core, rep(ncores, nlayer%%ncores))
    w <- lapply(1:ncores, function(x) list(i = x, data = w[[c(which(core == 
                                                                      x), which(core == x) + nlayer)]]))
    require(doParallel)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    wr <- foreach(x = w, .packages = c("raster", "windscape")) %dopar% 
      {
        x$data <- add_res(x$data)
        x$data <- add_lat(x$data)
        raster::calc(x$data, fun = rosefun, forceapply = TRUE, 
                     filename = paste0(substr(outfile, 1, nchar(outfile) - 
                                                4), "_", x$i, substr(outfile, nchar(outfile) - 
                                                                       3, nchar(outfile))))
      }
    stopCluster(cl)
  }
}




# generate downwind and upwind dispersal surfaces for one focal site 
origin <- matrix(c(-105, 40), ncol=2)
downwind <- conductance %>% 
  transition_stack(windflow, directions=8, symm=F, direction="downwind") %>%
  accCost(origin)
upwind <- conductance %>% 
  transition_stack(windflow, directions=8, symm=F, direction="upwind") %>%
  accCost(origin)

# convert seconds to hours, resturcture data, and plot
d <- stack(downwind / 3600, upwind / 3600) %>%
  rasterToPoints() %>%
  as.data.frame() %>%
  rename(downwind=layer.1, upwind=layer.2) %>%
  gather(direction, wind_hours, downwind, upwind)
ggplot() +
  facet_wrap(~direction) +
  geom_raster(data=d, aes(x, y, fill=wind_hours)) +
  geom_point(data=as.data.frame(origin), aes(V1, V2)) +
  scale_fill_gradientn(colors=c("yellow", "red", "blue", "black")) +
  theme_void() +
  theme(legend.position="top",
        strip.text=element_text(size=15))

