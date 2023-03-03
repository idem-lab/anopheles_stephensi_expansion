# inter-patch population spread modelling (dynamic range model on a network) toy example

# make fake data to play with
set.seed(2023-03-03)

# fake locations and populations of cities
n_locs <- 20
coords <- 100 * matrix(runif(n_locs * 2), ncol = 2)
pop <- rpois(n_locs, 10000 * rlnorm(n_locs, 3, 1))

# fake other covariates in these cities, and include the log human population
n_covs <- 3
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

# make fake time to invasion as a diffusion kernel form one of the cities
origin <- sample.int(n_locs, 1)
time_to_invasion <- ceiling(euclidean_distance[origin, ] / 2) + 1
n_times <- max(time_to_invasion)

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
diffusion_range <- normal(1, 1, truncation = c(0, Inf))
log_rel_diffusion_rate <- euclidean_distance / diffusion_range

# parametric gravity model dispersal rate on Euclidean distance with a scaling
# factor and power parameters on each component:
#   R_{o,d} = scaling * (P_o^a * P_d^b) / D_{o,d}^c
# where P_o is the population at origin, P_d is the population at destination,
# and a, b, and c are positive real-valued parameters. Without the scaling
# parameter and on the log scale this is:
#   a * log(P_o) + b * log(P_d) - c * log(D_{o,d})

pop_power_orig <- normal(1, 1, truncation = c(0, Inf))
pop_power_dest <- normal(1, 1, truncation = c(0, Inf))
dist_power <- normal(1, 1, truncation = c(0, Inf))

weight_orig <- log(pop) * pop_power_orig
weight_dest <- log(pop) * pop_power_dest
pop_weights <- weight_orig %*% t(weight_dest)

# add a 1 to the euclidean distance to fudge the diagonals to be non-zero
weight_dist <- log1p(euclidean_distance) * dist_power
log_rel_gravity_rate <- pop_weights - weight_dist

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
dispersal_intercept <- normal(0, 1)
diffusion_weight <- uniform(0, 1)
gravity_weight <- 1 - diffusion_weight
log_dispersal_matrix <- dispersal_intercept + gravity_weight * log_rel_gravity_rate + diffusion_weight * log_rel_diffusion_rate
dispersal_matrix <- exp(log_dispersal_matrix)

# now model suitability for the species as a site-specific growth rate, modeled
# as log-linear with all interactions

# create a matrix with the intercept and all interactions, then matrix multiply
# to model the location suitability (consider using a regularising prior on these to
# limit overfitting) - this will be used to model growth rates
cov_matrix <- model.matrix(~1 + a*b*c, covs)
n_coefs <- ncol(cov_matrix)
coefs <- normal(0, 1, dim = n_coefs)
loc_suitability <- cov_matrix %*% coefs
R <- exp(loc_suitability)
# now to model the dispersal and growth of populations over time

# Model *carrying capacity* of each location via covariates and population. We
# would expect this to be proportional to the population, given the same
# suitability. Use a Beverton-Holt type model here as in DSDM paper draft.

# Scale the mosquito population carrying capacity by the size of the city (a
# priori more people = more habitats, but the suitability model can adjust this
# using the covariate). Do this by defining a density-dependence effect
# parameter M which relates the growth rate to the carrying capacity in that city
alpha <- normal(0, 1)
M <- exp(alpha + log(pop))

# compute the carrying capacity - first raw according to Beverton Holt (negative
# if R < 1), and then rectified ~~0 if R < 1
K_raw <- expm1(loc_suitability) * M
K <- log1pe(K_raw)

# # # no longer doing this:
# # We need to constrain the sizes of the mosquito populations so that they are
# # not unidentified with other parameters. Cap one of the locations to be a
# # reference city with a fixed relative carrying capacity of 1, and the rest relative
# # to it
# reference_loc <- 1
# # K = rel_K / rel_K[reference_loc]
# log_K <- log_rel_K - log_rel_K[reference_loc]
# K <- exp(log_K)

# set an initial vector population state in these locations, at or near carrying
# capacity in native range
in_native_range <- seq_len(n_locs) == origin
initial_state <- in_native_range * K


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

transition_function <- function(state, iter, R, M, dispersal_matrix) {
  dd_scaling <- 1 + state / M
  grown_state <- (R * state) / dd_scaling
  dispersed_state <- dispersal_matrix %*% grown_state
  dispersed_state
}

# library(greta.dynamics)

# iters <- iterate_dynamic_function(transition_function = transition_function,
#                                 initial_state = initial_pop,
#                                 niter = n_times,
#                                 tol = 0,
#                                 R = R,
#                                 M = M,
#                                 dispersal_matrix = dispersal_matrix)

# # alternatively:
# matrix_function <- function(state, iter, R, M, dispersal_matrix) {
#   dd_scaling <- 1 + state / M
#   growth_rate <- R / dd_scaling
#   # guessing at the right dimension to sweep over here
#   new_matrix <- sweep(dispersal_matrix, 2, growth_rate, FUN = "*")
#   new_matrix
# }
# iters <- iterate_dynamic_matrix(matrix_function = matrix_function,
#                                 initial_state = initial_pop,
#                                 niter = n_times,
#                                 tol = 0,
#                                 R = R,
#                                 M = M,
#                                 dispersal_matrix = dispersal_matrix)

# states <- t(iters$all_states)

# solve with a for loop for now (slooow) until greta.dynamics is fixed up to
# work with TF2
states <- zeros(n_times, n_locs)
states[1, ] <- state <- initial_state

for (i in 2:n_times) {
  state <- transition_function(state, 1, R, M, dispersal_matrix)
  states[i, ] <- state
}

# slooooooow
calculate(states, nsim = 1)


# define the likelihood based on this the probability of detecting the species
# is modelled as the probability of detecting one or more individuals, under a
# poisson assumption of the number of individuals. Compute the likelihood of
# detection if the location had the highest level of surveillance, and then
# scale down from that with the observation probabilities (relative probability
# of detection). this doesn't need a scaling factor for the detection rate under
# perfect detection as that is confounded with the carrying capacity, so it's
# captured in the parameter alpha.

p_detect_perfect <- 1 - exp(-states)
p_detect <- p_detect_perfect * obs_probs

# define the likelihood
distribution(pa_data) <- bernoulli(p_detect)

# now we just need to fit it!

# 1. need to get greta.dynamics working on TF2

# 2. need to get initialisation hack working on TF2 (this is exponential growth
# so can go postal with numeric under/overflow)

# 3. need to constrain dispersal matrix to have all rows equal to some dispersal
# rate < 1, for identifiability between those parameters and alpha

# note: 
# this might cause the populations to grow from the initial very small
# dispersal, even in the furthest places. That might be fine, it's a limitation
# of the continuous model.




