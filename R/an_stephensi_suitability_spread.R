# Model An. stephensi spread and suitability via a dynamic range model with
# occupancy-detection likelihood. This is adapted from toy_spread_model.R to the
# data we have for An. stephensi

set.seed(2023-09-04)
library(tidyverse)
library(terra)
library(sf)
library(greta)
library(dismo)

# load data

# load pixel covariates
mask <- rast("output/rasters/derived/mask.tif")
larval_covs <- rast("output/rasters/derived/larval_habitat_covariates.tif")
populated <- rast("output/rasters/derived/populated_areas.tif")
climatic_rel_abund <- rast("output/rasters/derived/climatic_rel_abund.tif")

# An. stephensi presence, absence, and years of detection
first_detection <- readRDS("output/tabular/first_detection.RDS")

# filter this to points within the study region, and within the climatic
# suitability prediction
vals <- terra::extract(climatic_rel_abund, first_detection[, 1:2])[, 2]
first_detection <- first_detection %>%
  filter(is.finite(vals) & vals > 0)

# augment this with background points
non_detection <- dismo::randomPoints(
  mask = raster::raster(climatic_rel_abund),
  n = 500,
  p = as.matrix(first_detection[, 1:2])) %>%
  as_tibble() %>%
  mutate(
    year_first_detected = max(first_detection$year_first_detected),
    ever_detected = 0,
    background = 1
  )

# combine these, noting which are background points (lower detection
# probability), rather than absence records
detections <- bind_rows(
  non_detection,
  mutate(first_detection,
         background = 0)
)
# 
# plot(populated)
# points(detections,
#        pch = 21,
#        cex = 1,
#        bg = 1 - detections$background)

# load hexes
hexes <- readRDS("output/hexes.RDS")



# aggregate the values of these covariates at the hex scale


# note: we only need total K over a hex to model spread, so we could aggregate
# from pixel level to hex level at each MCMC iteration (logsumexp(log_K) over
# pixels in each hex), but for a prototype, just do the ecological fallacy thing
# and calculate summary stats for covariates over all urban-ish areas within a hex

log_climatic_suitability <- ?
larval_habitat_covs <- ?



plot(landcover)





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
larval_habitat_formula <- ~ -1 + built_env + pop + aridity
larval_habitat_cov_matrix <- model.matrix(larval_habitat_formula,
                                          data = larval_habitat_covs)
n_coefs <- ncol(larval_habitat_cov_matrix)
coefs <- normal(0, 1, dim = n_coefs)
log_larval_habitat_density <- cov_matrix %*% coefs

# get the carrying capacity
log_K <- log_climatic_suitability + log_larval_habitat_density
K <- exp(log_K)

# define the density-dependence scaling parameter M, expected to be positive,
# so a standard lognormal; with a prior median of 1
log_M <- normal(0, 1)
M <- exp(log_M)

# compute the log-intrinsic growth rate from this
log_R <- log1pe(log_K - log_M)

# set an initial vector population state in these locations, at or near carrying
# capacity in native range, and at some minimal level in other areas (note that
# masking with 0s broke the model gradients, for reasons I don't fully
# understand)
eps <- sqrt(.Machine$double.eps)
native_range_mask <- ifelse(df$in_native_range, 1, eps)
initial_state <- K * native_range_mask


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

transition_function <- function(state, iter, log_R, M, dispersal_matrix) {
  # dd_scaling <- 1 + state / M
  log_dd_scaling <- log1p(state / M)
  # fecundity <- R * state
  log_fecundity <- log_R + log(state)
  # grown_state <- fecundity / dd_scaling
  log_grown_state <- log_fecundity - log_dd_scaling
  grown_state <- exp(log_grown_state)
  dispersed_state <- dispersal_matrix %*% grown_state
  dispersed_state
}

library(greta.dynamics)

iters <- iterate_dynamic_function(transition_function = transition_function,
                                  initial_state = initial_state,
                                  niter = n_times,
                                  tol = 0,
                                  log_R = log_R,
                                  M = M,
                                  dispersal_matrix = dispersal_matrix)

# get relative abundance in hexes
rel_abund_hex <- t(iters$all_states)
log_rel_abund_hex <- log(rel_abund_hex)

# convert these to log proportions of carrying capacities in those hexes
log_prop_full_hex <- log_rel_abund_hex - log_K_hex

# multiply these by the pixel carrying capacities to get K there

# hexes to points is an integer lookup giving the hexes to which each of the
# points belong
log_rel_abund_points <- log_K_points + log_prop_full_hex[hexes_to_points]

# now define a likelihood for the (right truncated) time to detection dataset

# compute the relative probability of detection in a single year in a given
# location (relative to a density of 1 mosquito)
detection_probability <- normal(0, 0.1, truncation = c(0, 1))
# compute the probability of detection in that place and year
p_detect_annual <- 1 - exp(-rel_abund * detection_probability)
# and the cumulative probability of detection
p_detect_cumul <- apply(p_detect_annual, 1, "cumprob")

# model the right-censored time to detection with a time-dependent cox
# proportional hazards likelihood
log_detection_probability <- log(detection_probability)
log_hazard_annual <- log_rel_abund + log_detection_probability
log_hazard_cumul <- apply(log_hazard_annual, 1, "cumsum")
p_detect_cumul <- icloglog(log_hazard_cumul)

# define the Bernoulli likelihood over these times to detection
p_detect_annual <- 1 - exp(-exp(log_rel_abund + log_detection_probability))





# change this to have occupancy calculated on the point level (augmented with
# rangom background/integration points) - with the invasion process computed on
# hexes. I.e. occupancy of suitable pixels depends only on the hexes to which
# they belong.

# currently the invasion process depends on K (suitability), which is calculated
# as aggregated over all pixels.

# use integration points to sum K over cells?






