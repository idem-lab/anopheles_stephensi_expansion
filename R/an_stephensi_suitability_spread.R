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
populated <- rast("output/rasters/derived/populated_areas.tif")
climatic_rel_abund <- rast("output/rasters/derived/climatic_rel_abund.tif")

# add an epsilon to the climatic relative abundance to avoid 0 abundances later
climatic_rel_abund <- climatic_rel_abund + .Machine$double.eps
climatic_rel_abund <- mask(climatic_rel_abund, mask)

# An. stephensi presence, absence, and years of detection
first_detection <- readRDS("output/tabular/first_detection.RDS")

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

climatic_rel_abund_points <- terra::extract(climatic_rel_abund,
                                     detections[, 1:2],
                                     ID = FALSE) %>%
  as_tibble() %>%
  pull(mean)


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
larval_hab_formula <- ~ 1 + built_volume + pop

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
log_K_hex <- log(climatic_rel_abund_hexes) + log_larval_hab_hex
log_K_points <- log(climatic_rel_abund_points) + log_larval_hab_points
K_hex <- exp(log_K_hex)
K_points <- exp(log_K_points)

# define the density-dependence scaling parameter M, expected to be positive,
# so a standard lognormal; with a prior median of 1
log_M <- normal(0, 1)
M <- exp(log_M)

# compute the log-intrinsic growth rate from this
log_R_hex <- log1pe(log_K_hex - log_M)

# set an initial vector population state in these locations, at or near carrying
# capacity in native range, and at some minimal level in other areas (note that
# masking with 0s broke the model gradients, for reasons I don't fully
# understand)
# for all hexes, find out whether they contain at least one 'native range' point
native_range_lookup <- detections %>%
  group_by(id) %>%
  summarise(
    native = any(native)
  ) %>%
  right_join(
    tibble(
      id = hex_order),
    by = "id"
  ) %>%
  arrange(
    id
  ) %>%
  mutate(
    native = replace_na(native, FALSE)
  )
native_range_lookup <- native_range_lookup[match(native_range_lookup$id, hex_order),]
in_native_range <- native_range_lookup$native
# eps <- sqrt(.Machine$double.eps)
# native_range_mask <- ifelse(in_native_range, 1, eps)
# initial_state <- K_hex * native_range_mask
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

# compute the relative probability of detection in a single year in a given
# location (relative to a density of 1 mosquito)
detection_probability <- normal(0, 0.1, truncation = c(0, 1))
log_detection_probability <- log(detection_probability)
# compute the probability of detection in that place and year
p_detect_annual <- 1 - exp(-exp(log_rel_abund_points + log_detection_probability))
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
           detection_probability,
           log_M)

source("R/generate_valid_inits.R")
library(tensorflow)
n_chains <- 4
inits <- generate_valid_inits(m, n_chains)
draws <- mcmc(m,
              # warmup = 100,
              # n_samples = 100,
              chains = n_chains,
              initial_values = inits)

coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)
plot(draws)
summary(draws)

# make prediction rasters
all_cells <- seq_len(ncell(mask))
non_na_cells <- all_cells[!is.na(extract(mask, all_cells))]

# predict the larval habitat and carrying capacities
climatic_rel_abund_cells <- extract(climatic_rel_abund,
                                    non_na_cells)[, 1]
larval_covs_cells <- extract(larval_covs, non_na_cells) %>%
  as_tibble()
larval_hab_cov_cells <- model.matrix(larval_hab_formula,
                                     data = larval_covs_cells)
log_larval_hab_cells <- larval_hab_cov_cells %*% coefs
log_K_cells <- log(climatic_rel_abund_cells) + log_larval_hab_cells

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
                         # values = draws,
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

# quick plot of them
par(mfrow = c(3, 1))
plot(larval_habitat)
plot(K)
plot(max_occupancy)


# get posterior simulations
sims_temporal <- calculate(occupancy_cells,
                           # values = draws,
                           nsim = 100)

occupancy_cells_pred <- apply(sims_temporal$occupancy_cells,
                              2:3,
                              mean)

# set up time-varying rasters
mask_multi <- replicate(n_years, mask, simplify = FALSE) %>%
  do.call(c, .)
occupancy <- mask_multi
names(occupancy) <- years

for (i in seq_len(n_years)) {
  occupancy[[i]][non_na_cells] <- occupancy_cells_pred[i, ]
}

years_plot <- as.character(round(seq(min(years), max(years), length.out = 6)))
par(mfrow = c(3, 2))
plot(occupancy[[years_plot]])


# to do:

# add an observation probability parameter, based on the density of background
# points. Make this fixed, for now, and assume proportionality to the
# probability of detection

# think about adding in post-detection increases in observation probability for
# the same country and/or neighbouring countries

