# inter-patch population spread modelling (dynamic range model on a network) toy example

# make fake data to play with
set.seed(2023-03-03)

# fake locations and populations of cities
n_locs <- 20
coords <- 100 * matrix(runif(n_locs * 2), ncol = 2)
pop <- rpois(n_locs, 10000 * rlnorm(n_locs, 3, 1))

# fake other covariates in these cities, and include the log human population
n_covs <- 2
covs <- matrix(rnorm(n_locs * n_covs), n_locs, n_covs)
colnames(covs) <- letters[1:n_covs]
covs <- as.data.frame(cbind(covs, log(pop)))

# plot these
plot(coords,
     cex = 2 * pop / max(pop) + 1,
     pch = 16,
     col = scales::alpha("grey10", 0.5),
     asp = 1)

# build distance matrix
euclidean_distance <- as.matrix(dist(coords))

# make fake time to invasion as a diffusion kernel from one of the cities
n_times <- 10
origin <- sample.int(n_locs, 1)
time_to_invasion <- 1 + ceiling((n_times - 1) * euclidean_distance[origin, ] / max(euclidean_distance[origin, ]))
# n_times <- max(time_to_invasion)

# plot these
pal <- colorRampPalette(c("blue", "grey90"))
plot(coords,
     cex = 2 * pop / max(pop) + 1,
     pch = 16,
     col = pal(n_times)[time_to_invasion],
     asp = 1)

# make a 'true' presence/absence dataset, n_times x n_locs
true_presence <- sapply(time_to_invasion,
                        function(time) {
                          (1:n_times) >= time
                        })

# add observation probabilities for each site, biased by population and
# increasing over time with awareness
site_logit_obs_prob <- (log(pop) - log(mean(pop))) / 5
time_logit_obs_prob <- log(2 + 1:n_times)
obs_probs <- plogis(outer(time_logit_obs_prob, site_logit_obs_prob, FUN = "+"))

# simulate some fake detection data
pa_data <- matrix(
  rbinom(n_times * n_locs,
         1,
         c(true_presence * obs_probs)),
  n_times, n_locs)

# Fit greta model with a modelled mobility matrix, city-level growth rate model,
# iterated with greta.dynamics

# No repeat visits in the same season, so we can't identify the detection
# probability, instead we need to assume it is fixed (and correct), so use the
# 'true' value above, incorporated with the probability of presence.

library(greta)
# build a mobility matrix model from two components:
# euclidean distance with a dispersal kernel, a gravity model.

# Exponential diffusion kernel on euclidean distance. The rate of movement
# (R_{o,d}) between the origin (o) and destination (d) is given by:
#   R_{o,d} = scaling * exp(D_{o,d} / diffusion_range)
# where D_{o,d} is the Euclidean distance between o and d, and dispersal_range
# is a positive real-valued model parameter controlling the average dispersal
# distance. We define this as scale-less (since that is confounded with the
# gravity model scale and we add it later) and on the log scale:

# average dispersal range
diffusion_range <- normal(1, 1, truncation = c(0, Inf))
# hist(calculate(diffusion_range, nsim = 1000)[[1]], breaks = 100)

# log relative rate of dispersal
log_unscaled_diffusion_rate <- euclidean_distance / diffusion_range

# relative rate of dispersal (including zero-length dispersal)
unscaled_diffusion_matrix <- exp(log_unscaled_diffusion_rate)

# set the diagonal elements to 0 (don't model rate of same-patch dispersal with
# this exponential)
diagonal <- diag(nrow = n_locs, ncol = n_locs)
unscaled_diffusion_matrix <- unscaled_diffusion_matrix * (1 - diagonal)

# make these columns sum to 1 to get probability of moving to each other patch *if* they left
rel_dispersal_matrix <- sweep(unscaled_diffusion_matrix, 1, rowSums(unscaled_diffusion_matrix), FUN = "/")

# probability of dispersion
dispersal_fraction <- normal(0, 0.5, truncation = c(0, 1))
# range(calculate(dispersal_fraction, nsim = 10000)[[1]])

# normalise these to have the overall probability of dispersing to that patch,
# and add back the probability of remaining
dispersal_matrix <- dispersal_fraction * rel_dispersal_matrix + (1 - dispersal_fraction) * diagonal

# # check this doesn't lead to population growth
# rowSums(calculate(dispersal_matrix, nsim = 1)[[1]][1, , ])
# growth <- dispersal_matrix %*% ones(n_locs)
# calculate(growth, nsim = 1)[[1]][1, , ]

# parametric gravity model dispersal rate on Euclidean distance with a scaling
# factor and power parameters on each component:
#   R_{o,d} = scaling * (P_o^a * P_d^b) / D_{o,d}^c
# where P_o is the population at origin, P_d is the population at destination,
# and a, b, and c are positive real-valued parameters. Without the scaling
# parameter and on the log scale this is:
#   a * log(P_o) + b * log(P_d) - c * log(D_{o,d})

# pop_power_orig <- normal(1, 1, truncation = c(0, Inf))
# pop_power_dest <- normal(1, 1, truncation = c(0, Inf))
# dist_power <- normal(1, 1, truncation = c(0, Inf))
# 
# weight_orig <- log(pop) * pop_power_orig
# weight_dest <- log(pop) * pop_power_dest
# pop_weights <- weight_orig %*% t(weight_dest)
# 
# # add a 1 to the euclidean distance to fudge the diagonals to be non-zero
# weight_dist <- log1p(euclidean_distance) * dist_power
# log_rel_gravity_rate <- pop_weights - weight_dist

# combine these into a dispersal matrix. We could apply a separate scaling to
# each one, sum them, and then exponentiate:
#   gravity_weight <- normal(0, 1)
#   diffusion_weight <- normal(0, 1)
#   log_dispersal_matrix <- gravity_weight * log_rel_gravity_rate + diffusion_weight * log_rel_diffusion_rate
#   dispersal_matrix <- exp(log_dispersal_matrix)

# however since these log_rel_*_rate matrices will be correlated (distance is a
# major driver and appears in both), the matrix weights will be correlated,
# which will lead to bad sampling. Instead, we can decompose them into an
# overall weight and a proportion assigned to the components.

# log_dispersal_matrix <- log(dispersal_fraction) * log_rel_diffusion_rate
# gravity_weight <- 1 - diffusion_weight
# log_dispersal_matrix <- dispersal_intercept + gravity_weight * log_rel_gravity_rate + diffusion_weight * log_rel_diffusion_rate


# now model suitability for the species as a site-specific growth rate, modeled
# as log-linear with all interactions

# create a matrix with the intercept and all interactions, then matrix multiply
# to model the location suitability (consider using a regularising prior on these to
# limit overfitting) - this will be used to model growth rates
cov_matrix <- model.matrix(~1 + a*b, covs)
n_coefs <- ncol(cov_matrix)
# suitability_intercept <- normal(1, 1)
coefs <- normal(0, 1, dim = n_coefs)
log_R <- cov_matrix %*% coefs
R <- exp(log_R)
# now to model the dispersal and growth of populations over time

# Model *carrying capacity* of each location via covariates and population. We
# would expect this to be proportional to the population, given the same
# suitability. Use a Beverton-Holt type model here as in DSDM paper draft.

# Scale the mosquito population carrying capacity by the size of the city (a
# priori more people = more habitats, but the suitability model can adjust this
# using the covariate). Do this by defining a density-dependence effect
# parameter M which relates the growth rate to the carrying capacity in that city
alpha <- normal(-10, 1)
M <- exp(alpha + log(pop))

# compute the carrying capacity - first raw according to Beverton Holt (negative
# if R < 1), and then rectified ~~0 if R < 1
K_raw <- expm1(log_R) * M
K <- log1pe(K_raw)

# calculate(R, nsim = 1)[[1]][1, , ]
# tmp <- calculate(K, nsim = 1)[[1]][1, , ]
# round(tmp)

# set an initial vector population state in these locations, at or near carrying
# capacity in native range, and at some minimal level in other areas (note that
# masking with 0s broke the model gradients, for reasons I don't fully
# understand)
in_native_range <- seq_len(n_locs) == origin
eps <- sqrt(.Machine$double.eps)
native_range_mask <- ifelse(in_native_range, 1, eps)
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
# bit) implemented in log space for a) numerical stability and b) because fo a
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
states <- t(iters$all_states)

# tmp <- calculate(states, nsim = 1)[[1]][1, , ]
# round(tmp, 0)

# # solve with a for loop for now (slooow) until greta.dynamics is fixed up to
# # work with TF2
# states <- zeros(n_times, n_locs)
# states[1, ] <- state <- initial_state
# 
# for (i in 2:n_times) {
#   state <- transition_function(state, 1, R, M, dispersal_matrix)
#   states[i, ] <- state
# }
# # 
# # # slooooooow
# # calculate(states, nsim = 1)


# define the likelihood based on this the probability of detecting the species
# is modelled as the probability of detecting one or more individuals, under a
# poisson assumption of the number of individuals. Compute the likelihood of
# detection if the location had the highest level of surveillance, and then
# scale down from that with the observation probabilities (relative probability
# of detection). this doesn't need a scaling factor for the detection rate under
# perfect detection as that is confounded with the carrying capacity, so it's
# captured in the parameter alpha.

# mock up models without more complex bits
# states <- sweep(ones(n_times, n_locs), 2, initial_state, FUN = "*")

p_detect_perfect <- 1 - exp(-states)
p_detect <- p_detect_perfect * obs_probs

# tmp <- calculate(p_detect[n_times, ], nsim = 1)[[1]][1, , ]
# round(tmp, 3)
# pa_data[n_times, ]

# define the likelihood
distribution(pa_data) <- bernoulli(p_detect)

m <- model(alpha,
           coefs,
           dispersal_fraction,
           diffusion_range)

# get good starting values for chains (difficult because of exponential growth)
n_chains <- 30
source("R/generate_valid_inits.R")
library(tensorflow)
inits <- generate_valid_inits(m, n_chains)
draws <- mcmc(m, chains = n_chains, initial_values = inits)

plot(draws)
coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)

bayesplot::mcmc_pairs(draws)
# plot(calculate(states[1, 1], values = draws))
# summary(calculate(states[1, 1], values = draws, nsim = 1000)[[1]])

# now need to work out identifiability for multiple component dispersal
# matrices:
# - model each one, and constrain it to have rows summing to 1, then define a
# set of weights (summing to 1) across the matrices
# - need to make sure there is no intercept parameter in each that will be
# non-identified with the weights

# make predictions of time to invasion
present_sim <- calculate(p_detect_perfect > 0.5,
                         values = draws,
                         nsim = 1000)[[1]]

first <- function(x) min(c(which(as.logical(x)), length(x)))
arrival_time_sim <- apply(present_sim, c(1, 3), first)

mean_arrival_time <- colMeans(arrival_time_sim)
ci_arrival_time <- apply(arrival_time_sim, 2, quantile, c(0.025, 0.975))

# plot true and estimated arrival times
par(mfrow = c(2, 1))
plot(coords,
     cex = 2 * pop / max(pop) + 1,
     pch = 16,
     col = pal(n_times)[time_to_invasion],
     asp = 1)
points(x = coords[origin, 1], y = coords[origin, 2],
       cex = 2 * pop[origin] / max(pop) + 1,
       pch = 21,
       lwd = 2,
       col = "red")
title(main = "truth")

plot(coords,
     cex = 2 * pop / max(pop) + 1,
     pch = 16,
     col = pal(n_times)[ceiling(mean_arrival_time)],
     asp = 1)
points(x = coords[origin, 1], y = coords[origin, 2],
       cex = 2 * pop[origin] / max(pop) + 1,
       pch = 21,
       lwd = 2,
       col = "red")
title(main = "estimated")




# Note that the fitted model is different from the 'true' generating process for
# this simple model. The fitted model accounts for large suitable patches acting
# as sources in the spread process. So it's not bad considering that.
