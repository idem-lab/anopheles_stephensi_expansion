# Model An. stephensi spread and suitability via a dynamic range model with
# occupancy-detection likelihood. This is adapted from toy_spread_model.R to the
# data we have for An. stephensi

set.seed(2023-09-04)
library(tidyverse)
library(terra)
library(sf)
library(greta)
library(dismo)

# load all functions
. <- list.files(path = "R/functions/",
                pattern = "*.R",
                full.names = TRUE) %>%
  lapply(source)

# load data

# load hexes and hex raster lookup (all cells and populated areas)
hexes <- readRDS("output/hexes.RDS")
hex_lookup <- rast("output/rasters/derived/hex_raster_lookup.tif")
hex_populated_lookup <- rast("output/rasters/derived/hex_raster_populated_lookup.tif")

# load pixel covariates
mask <- rast("output/rasters/derived/mask.tif")
larval_covs <- rast("output/rasters/derived/larval_habitat_covariates.tif")
climatic_rel_abund <- rast("output/rasters/derived/climatic_rel_abund.tif")
as_detection_density <- rast("output/rasters/derived/an_stephensi_detection_density.tif")

mask <- mask(mask, hex_lookup)
hex_populated_lookup <- mask(hex_populated_lookup, mask)

# make a raster of cell areas, and scale them
cell_area <- terra::cellSize(mask)
cell_area <- mask(cell_area, mask)
cell_area <- cell_area / global(cell_area, "max", na.rm = TRUE)[1, 1]

# add an epsilon to the climatic relative abundance to avoid 0 abundances later
climatic_rel_abund <- mask(climatic_rel_abund, mask)
climate_suitable <- climatic_rel_abund > 0
climatic_rel_abund <- climatic_rel_abund + .Machine$double.eps

# An. stephensi presence, absence, and years of detection
first_detection <- readRDS("output/tabular/first_detection.RDS")

# move some of these points onto land (up to 50km), and remove the remainder
# (some outside the region, others on tiny islands)
first_detection <- first_detection %>%
  mutate(
    nearest_land(
      coords = dplyr::select(., x, y),
      raster = mask,
      projected_crs = "ESRI:102022",
      max_distance_m = 50 * 1e3
    )
  ) %>%
  filter(
    !is.na(x) & !is.na(y)
  )

# filter out any An stephensi detections that are inconsistent with the climatic
# relative abundance (parts of Iran, Afghanistan)
first_detection <- first_detection %>%
  mutate(
    valid = !is.na(terra::extract(mask, .[, c("x", "y")])[, 2]),
    suitable = terra::extract(climate_suitable, .[, c("x", "y")])[, 2],
    suitable = case_when(
      !is.finite(suitable) ~ FALSE,
      .default = suitable)
  )

# filter this to points within the study region, and within the climatic
# suitability prediction
climate_vals <- terra::extract(climatic_rel_abund, first_detection[, 1:2])[, 2]
hex_vals <- terra::extract(hex_populated_lookup, first_detection[, 1:2])[, 2]
first_detection <- first_detection %>%
  filter(is.finite(hex_vals) & climate_vals > 0)

# augment this with background points
non_detection <- dismo::randomPoints(
  mask = raster::raster(hex_populated_lookup),
  n = 500,
  p = as.matrix(first_detection[, 1:2])) %>%
  as_tibble() %>%
  mutate(
    native = FALSE,
    year_first_detected = max(first_detection$year_first_detected),
    ever_detected = 0,
    background = 1
  )

# combine these, noting which are background points (lower detection
# probability), rather than absence records
detections <- non_detection %>%
  bind_rows(
    mutate(first_detection,
           background = 0)
  ) %>%
  mutate(
    id = terra::extract(hex_populated_lookup,
                        .[, 1:2],
                        ID = FALSE)[, 1]
  )

# get all years
years <- seq(min(detections$year_first_detected),
             max(detections$year_first_detected))
n_years <- length(years)

# extract the values of the larval covariates at the points
larval_covs_points <- terra::extract(larval_covs,
                                     detections[, 1:2],
                                     ID = FALSE) %>%
  as_tibble()

# and the climatic model
climatic_rel_abund_points <- terra::extract(climatic_rel_abund,
                                     detections[, 1:2],
                                     ID = FALSE) %>%
  as_tibble() %>%
  pull(mean)

# and the detection covariates (log(1+x) transformed)
log_detection_density_points <- terra::extract(as_detection_density,
                                               detections[, 1:2],
                                               ID = FALSE) %>%
  as.matrix() %>%
  t() %>%
  log1p()

# and the log cell areas for points
log_area_points <- terra::extract(cell_area,
                                  detections[, 1:2],
                                  ID = FALSE) %>%
  as_tibble() %>%
  pull(area) %>%
  log()

# aggregate the values of the covariates at the hex scale
covs_continuous <- !names(larval_covs) %in% c("type", "cover")
larval_covs_continuous_hex <- terra::zonal(larval_covs[[covs_continuous]],
                                           hex_populated_lookup,
                                           fun = "mean",
                                           na.rm = TRUE) %>%
  as_tibble()

larval_covs_discrete_hex <- terra::zonal(larval_covs[[!covs_continuous]],
                                           hex_populated_lookup,
                                           fun = "modal",
                                           na.rm = TRUE) %>%
  as_tibble()

larval_covs_hex <- larval_covs_continuous_hex %>%
  left_join(larval_covs_discrete_hex,
            by = "id")
hex_order <- larval_covs_hex %>%
  pull(id)

larval_covs_hex <- larval_covs_hex %>%
  dplyr::select(-id)

# force the same order
larval_covs_hex <- larval_covs_hex[, names(larval_covs_points)]

# extract the climatic relative abundances too
climatic_rel_abund_hexes <- terra::zonal(
  x = climatic_rel_abund,
  z = hex_populated_lookup,
  fun = "mean") %>%
  as_tibble() %>%
  pull(mean)

# and the total cell areas (then divide by max to maintain scaling on M)
log_area_hex <- terra::zonal(
  x = cell_area,
  z = hex_populated_lookup,
  fun = "sum")  %>%
  as_tibble() %>%
  pull(area) %>%
  log()
# area_hexes <- area_hexes / max(area_hexes)
log_area_hex <- log_area_hex - max(log_area_hex)

# load a dispersal matrix between hexes
radiation_dispersal_raw <- readRDS("output/tabular/radiation_matrix.RDS")
gravity_dispersal_raw <- readRDS("output/tabular/gravity_matrix.RDS")
distance_matrix_raw <- readRDS("output/tabular/distance_matrix.RDS")

# create a lookup to the hex ids, to keep only those we use
matrix_order <- tibble(
  index = seq_len(nrow(radiation_dispersal_raw)),
  h3_index = rownames(radiation_dispersal_raw)
)

matrix_index <- hexes %>%
  left_join(
    matrix_order,
    by = "h3_index") %>% 
  pull(index)

n_hexes <- nrow(hexes)

# subset all of these
radiation_dispersal <- radiation_dispersal_raw[matrix_index, matrix_index]
gravity_dispersal <- gravity_dispersal_raw[matrix_index, matrix_index]
distance_matrix <- distance_matrix_raw[matrix_index, matrix_index]

# model exponential dispersal kernel with the distance matrix, (accounting for
# different ranges) and zero out the diagonals
distance_decay <- normal(0, 0.01, truncation = c(0, Inf))
exponential_dispersal <- exp(-distance_matrix / distance_decay)

# overall rate of dispersal
dispersal_fraction <- normal(0, 0.01, truncation = c(0, 1))

# apply weights (summing to 1) for the different dispersal modes, with a priori
# equal weights on each
# dispersal_weights <- t(dirichlet(alpha = t(rep(1/3, 3))))

# don't use dirichlet, because it's hard to initialise
multilogit_weights <- normal(0, 1, dim = c(1, 2))
dispersal_weights <- t(imultilogit(multilogit_weights))
# hist(calculate(dispersal_weights[2], nsim = 1000)[[1]], breaks = 100)

weighted_dispersal <- dispersal_weights[1] * radiation_dispersal +
  dispersal_weights[2] * gravity_dispersal +
  dispersal_weights[3] * exponential_dispersal

# apply the dispersal fraction to get the importation rate into each cell (not
# removing from the origin cell)
dispersal_matrix <- dispersal_fraction * weighted_dispersal

# set the diagonal elements to 1 so the matrix multiply doesn't affect
# within-location population
diag(dispersal_matrix) <- ones(n_hexes)

# model the observation probability with an intercept and positive-constrained
# coefficient for density, on the logit scale
det_prob_intercept <- normal(qlogis(0.01), 1)
det_prob_slope <- normal(0, 1, truncation = c(0, Inf))

# model the relative probability of detection in a single year in a given
# location (relative to a density of 1 mosquito)
logit_det_prob <- det_prob_intercept + det_prob_slope * log_detection_density_points
# I *think* this is more numerically stable:
log_det_prob <- log(ilogit(logit_det_prob))
# log_det_prob <- -log1p(exp(-logit_det_prob))
# log(ilogit(x))
# = log(1 / (1 + exp(-x)))
# = log(1) - log(1 + exp(-x))
# = -log1p(exp(-x))

# build model

# note: we only need total K over a hex to model spread, so we could aggregate
# from pixel level to hex level at each MCMC iteration (logsumexp(log_K) over
# pixels in each hex), but for a prototype, just do the ecological fallacy thing
# and calculate summary stats for covariates over all urban-ish areas within a
# hex, and model with those

# build up a pixel-level model for the carrying capacity for An. stephensi

# from Beverton-Holt, the carrying-capacity K is:
#   K = (R - 1) * M
# where R is the intrinsic growth rate of the population, and M is a constant
# representing the density-dependent effect. We model K directly from adult
# abundance relative to larval habitat, and larval habitat density, and can
# infer M, so we rearrange to calculate R as:
#   R = 1 + K / M
# on the log scale we have:
#   log(R) = log(1 + K / M)
#   log(R) = log1p(K / M)
#   log(R) = log1pe(log(K) - log(M))

# model the relative density of larval habitats (note the intercept term here
# controls the absolute value of K, so remove it to make it relative and control
# that instead in the observation probabilities)
larval_hab_formula <- ~ -1 + built_volume + tcb + tcw

# extract for both points and hexes
larval_hab_cov_points <- model.matrix(larval_hab_formula,
                                      data = larval_covs_points)
larval_hab_cov_hex <- model.matrix(larval_hab_formula,
                                      data = larval_covs_hex)
n_coefs <- ncol(larval_hab_cov_points)
coefs <- normal(0, 1, dim = n_coefs)

# get log densities for points and hexes
log_larval_hab_points <- larval_hab_cov_points %*% coefs
log_larval_hab_hex <- larval_hab_cov_hex %*% coefs

# get the carrying capacity at hex and point scale
log_K_hex <- log(climatic_rel_abund_hexes) +
  log_larval_hab_hex +
  log_area_hex
log_K_points <- log(climatic_rel_abund_points) +
  log_larval_hab_points +
  log_area_points
K_hex <- exp(log_K_hex)
K_points <- exp(log_K_points)

# define the density-dependence scaling parameter M, expected to be positive,
# so a standard lognormal; with a prior median of 1
log_M <- normal(0, 1)
M <- exp(log_M)

# compute the log-intrinsic growth rate from this
log_R_hex <- log1pe(log_K_hex - log_M)

# manually define the hexes in the native range to roughly match that in Sinka
# et al 2020 (exclude western part of Arabian peninsula)
hex_native <- c(7, 8, 11, 27, 29, 38, 39, 40, 41, 45, 46, 54, 
              68, 69, 71, 72, 73, 74, 75,
              93, 94, 99, 100, 101, 102)
in_native_range <- as.numeric(hex_order %in% hex_native)

eps <- sqrt(.Machine$double.eps)
initial_state <- K_hex * in_native_range + eps * (1 - in_native_range)

# define the function to iterate dynamics

# the absolute carrying capacity is relative to detectability! So model
# reproduction numbers instead and e.g. Beverton-Holt with softplus and an
# intercept to map to the relative carrying capacity. No need to normalise.

# update populations as:
#   N_{t+1} = (R * N_{t}) / (1 + N_{t} / M)

# where M is a positive real density-dependence effect which we can model as a
# single scaling factor / intercept term. We get K from M and R as:
#   K = (R - 1) *  M
#   M = K / (R - 1)

# at each iteration:
# - update the population in each location based on the carrying capacity (force a link between growth rate and carrying capacity)
# - then do dispersal

# this is a Beverton-Holt population dynamic model (density dependent growth
# bit) implemented in log space for a) numerical stability and b) because of a
# bug in greta.dynamics that makes it fail with the '1' specified inside the
# function
transition_function <- function(state, iter, K, log_R, M, dispersal_matrix) {
  # dd_scaling <- 1 + state / M
  log_dd_scaling <- log1p(state / M)
  # fecundity <- R * state
  log_fecundity <- log_R + log(state)
  # grown_state <- fecundity / dd_scaling
  log_grown_state <- log_fecundity - log_dd_scaling
  grown_state <- exp(log_grown_state)
  # now compute importations into each cell, from the dispersal matrix (with
  # diagonal elements set to 1 so no within-cell population change)
  dispersed_state <- dispersal_matrix %*% grown_state
  # cap this at the carrying capacity and return
  capped_dispersed_state <- greta_pmin(dispersed_state, K) 
  capped_dispersed_state
}

library(greta.dynamics)

iters <- iterate_dynamic_function(transition_function = transition_function,
                                  initial_state = initial_state,
                                  niter = n_years,
                                  tol = 0,
                                  K = K_hex,
                                  log_R = log_R_hex,
                                  M = M,
                                  dispersal_matrix = dispersal_matrix)

# get relative abundance in hexes
rel_abund_hex <- t(iters$all_states)
log_rel_abund_hex <- log(rel_abund_hex)

# convert these to log proportions of carrying capacities in those hexes
log_prop_full_hex <- sweep(log_rel_abund_hex, 2, log_K_hex, FUN = "-")

# multiply these by pixel carrying capacities to get the relative abundance in
# the pixels

# hexes to points is an integer lookup giving the hexes to which each of the
# points belong
log_prop_full_points <- log_prop_full_hex[, detections$id]
log_rel_abund_points <- sweep(log_prop_full_points, 2, log_K_points, FUN = "+") 

# now define a likelihood for the (right truncated) time to detection dataset

# model detection probability, relative to weighted cumulative number of
# detections in the nearby area

# compute the probability of detection in that place and year
p_detect_annual <- 1 - exp(-exp(log_rel_abund_points + log_det_prob))
# and the cumulative probability of detection
p_detect_cumul <- 1 - t(apply(1 - p_detect_annual, 1, "cumprod"))

# can we instead model the right-censored time to detection with a
# time-dependent Cox proportional hazards likelihood to make this cumprod a
# cumsum?

# extract the relevant element for the censored likelihood
year_index <- match(detections$year_first_detected, years)
obs_index <- cbind(year_index, seq_len(nrow(detections)))

obs_prob <- p_detect_cumul[obs_index]

distribution(detections$ever_detected) <- bernoulli(obs_prob)

m <- model(coefs,
           dispersal_fraction,
           dispersal_weights,
           distance_decay,
           det_prob_intercept,
           det_prob_slope,
           log_M)

# fit the model
source("R/generate_valid_inits.R")
library(tensorflow)
n_chains <- 10
inits <- generate_valid_inits(m, n_chains)
draws <- mcmc(m,
              chains = n_chains,
              initial_values = inits)

coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

plot(draws)
summary(draws)

# identify the chains that fell in each mode

mn <- lapply(draws, function(x) mean(x[, "log_M"]))
big_m_chain <- mn > 2
mode_1_index <- which(big_m_chain)
mode_2_index <- which(!big_m_chain)

# need to split the draws, and secret draws, by these
mi <- attr(draws, "model_info")
mi_mode_1 <- mi_mode_2 <- mi

mi_mode_1$raw_draws <- mi$raw_draws[mode_1_index]
draws_mode_1 <- greta:::as_greta_mcmc_list(
  x = draws[mode_1_index],
  model_info = mi_mode_1
)

mi_mode_2$raw_draws <- mi_mode_2$raw_draws[mode_2_index]
draws_mode_2 <- greta:::as_greta_mcmc_list(
  x = draws[mode_2_index],
  model_info = mi_mode_2
)


# make prediction rasters

all_cells <- seq_len(ncell(mask))
non_na_cells <- all_cells[!is.na(extract(mask, all_cells))]

# predict the larval habitat and carrying capacities
climatic_rel_abund_cells <- extract(climatic_rel_abund,
                                    non_na_cells)[, 1]
log_area_cells <- extract(cell_area,
                          non_na_cells)[, 1] %>%
  log()
larval_covs_cells <- extract(larval_covs, non_na_cells) %>%
  as_tibble()
larval_hab_cov_cells <- model.matrix(larval_hab_formula,
                                     data = larval_covs_cells)
log_larval_hab_cells <- larval_hab_cov_cells %*% coefs
log_K_cells <- log(climatic_rel_abund_cells) +
  log_larval_hab_cells +
  log_area_cells

larval_hab_cells <- exp(log_larval_hab_cells)
K_cells <- exp(log_K_cells)

# and maximum occupancy
max_occupancy_cells <- 1 - exp(-K_cells)

# get the log proportions full in each timestep (hex level) and use it to
# compute relative abundance at cell level
log_prop_full_cells <- log_prop_full_hex[, hex_lookup[non_na_cells][, 1]]
log_rel_abund_cells <- sweep(log_prop_full_cells, 2, log_K_cells, FUN = "+") 
rel_abund_cells <- exp(log_rel_abund_cells)

# convert to occupancy
occupancy_cells <- 1 - exp(-rel_abund_cells)

# get posterior simulations for temporally static variables
sims_static <- calculate(K_cells,
                         larval_hab_cells,
                         max_occupancy_cells,
                         values = draws,
                         nsim = 100)

# make static maps
K_cells_pred <- apply(sims_static$K_cells, 2, mean)
larval_hab_cells_pred <- apply(sims_static$larval_hab_cells, 2, mean)
max_occupancy_cells_pred <- apply(sims_static$max_occupancy_cells, 2, mean)

larval_habitat <- K <- max_occupancy <- mask
names(K) <- "carrying capacity"
names(max_occupancy) <- "potential distribution"
names(larval_habitat) <- "larval habitat"

K[non_na_cells] <- K_cells_pred
larval_habitat[non_na_cells] <- larval_hab_cells_pred
max_occupancy[non_na_cells] <- max_occupancy_cells_pred


# make a binary layer of populated areas, to set all the non-populated areas to
# 0 in predictions (populated_areas.tif has the wrong extent - need to fix and align
# all of them!)
populated_area_mask <- hex_populated_lookup * 0 + 1
populated_area_mask[is.na(populated_area_mask)] <- 0
populated_area_mask <- mask(populated_area_mask, mask)

larval_habitat <- larval_habitat * populated_area_mask
K <- K * populated_area_mask
max_occupancy <- max_occupancy * populated_area_mask

# quick plot of them
par(mfrow = c(3, 1))
plot(larval_habitat)
plot(K)
plot(max_occupancy)

par(mfrow = c(1, 1))
plot(max_occupancy)

# save this out
terra::writeRaster(
  x = max_occupancy,
  filename = "output/rasters/derived/predicted_potential_distribution.tif",
  overwrite = TRUE
)

# get posterior simulations
sims_temporal <- calculate(occupancy_cells,
                                  values = draws,
                                  nsim = 100)
sims_temporal_mode_1 <- calculate(occupancy_cells,
                                  values = draws_mode_1,
                                  nsim = 100)
sims_temporal_mode_2 <- calculate(occupancy_cells,
                                  values = draws_mode_2,
                                  nsim = 100)

occupancy_cells_pred <- apply(sims_temporal$occupancy_cells,
                                     2:3,
                                     mean)
occupancy_cells_pred_mode_1 <- apply(sims_temporal_mode_1$occupancy_cells,
                                     2:3,
                                     mean)
occupancy_cells_pred_mode_2 <- apply(sims_temporal_mode_2$occupancy_cells,
                                     2:3,
                                     mean)

# set up time-varying rasters
mask_multi <- replicate(n_years, mask, simplify = FALSE) %>%
  do.call(c, .)
occupancy <- mask_multi
names(occupancy) <- years

occupancy_mode_1 <- occupancy_mode_2 <- occupancy

for (i in seq_len(n_years)) {
  occupancy[[i]][non_na_cells] <- occupancy_cells_pred[i, ]
  occupancy_mode_1[[i]][non_na_cells] <- occupancy_cells_pred_mode_1[i, ]
  occupancy_mode_2[[i]][non_na_cells] <- occupancy_cells_pred_mode_2[i, ]
}

# set unpopulated areas to 0
occupancy <- occupancy * populated_area_mask
occupancy_mode_1 <- occupancy_mode_1 * populated_area_mask
occupancy_mode_2 <- occupancy_mode_2 * populated_area_mask

years_plot <- as.character(round(seq(min(years), max(years), length.out = 6)))

years_plot <- c("1990", "2000", "2010", "2022")
par(mfrow = c(2, 2))
plot(occupancy[[years_plot]])

par(mfrow = c(2, 2))
plot(occupancy_mode_1[[years_plot]])

par(mfrow = c(2, 2))
plot(occupancy_mode_2[[years_plot]])

# save the distribution rasters
terra::writeRaster(
  x = occupancy,
  filename = "output/rasters/derived/predicted_occupancy.tif",
  overwrite = TRUE
)

terra::writeRaster(
  x = occupancy_mode_1,
  filename = "output/rasters/derived/predicted_occupancy_mode_1.tif",
  overwrite = TRUE
)

terra::writeRaster(
  x = occupancy_mode_2,
  filename = "output/rasters/derived/predicted_occupancy_mode_2.tif",
  overwrite = TRUE
)


# output proportion full, by full hexes, at least for debugging

prop_full_cells <- exp(log_prop_full_cells)

sims_temporal2 <- calculate(prop_full_cells,
                            values = draws_mode_2,
                            nsim = 100)

prop_full_cells_pred <- apply(sims_temporal2$prop_full_cells,
                              2:3,
                              mean)

# set up time-varying rasters
prop_full <- mask_multi
names(prop_full) <- years

for (i in seq_len(n_years)) {
  prop_full[[i]][non_na_cells] <- prop_full_cells_pred[i, ]
}

# set unpopulated areas to 0
# prop_full <- prop_full * populated_area_mask

years_plot <- as.character(round(seq(min(years), max(years), length.out = 6)))

par(mfrow = c(3, 2))
plot(prop_full[[years_plot]],
     range = c(0, 1))

# save the distribution rasters
terra::writeRaster(
  x = prop_full,
  filename = "output/rasters/derived/predicted_prop_full.tif",
  overwrite = TRUE
)
