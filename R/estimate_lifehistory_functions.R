# Estimate life history functions (life history parameters as a function of
# environmental conditions) from a variety of sources, including:
# Villena et al. https://doi.org/10.1002/ecy.3685

# First: refit the models in Villena et al. to get and interpolate the posterior
# predicted relationship against temperature of key life history parameters.
# This requires adapting the R and JAGS code provided in Villena et al. (the
# basis of the below code), and defining here a number of functions that code
# appears to use from an R package that was formerly hosted on github:
# lorecatta/DENVclimate, which appears to no longer exist anymore, but the
# source code for which is still available on rdrr.io

# These functions are defined below, sourced, reformatted, and slightly modified
# from: https://rdrr.io/github/lorecatta/DENVclimate/src/R/mcmc_utils_all.R
make.briere.samps <- function(coda.samps,
                              nchains = 2,
                              samp.lims = c(151, 5000),
                              sig = TRUE) {
  T0 <- Tm <- cc <- sigma <- NULL
  l1 <- samp.lims[1]
  l2 <- samp.lims[2]
  for (i in 1:nchains) {
    T0 <- c(T0, coda.samps[[i]][l1:l2, 1])
    Tm <- c(Tm, coda.samps[[i]][l1:l2, 2])
    cc <- c(cc, coda.samps[[i]][l1:l2, 3])
    if (sig) sigma <- c(sigma, coda.samps[[i]][l1:l2, 4])
  }
  if (sig) {
    samps <- data.frame(matrix(c(T0, Tm, cc, sigma),
                               ncol = 4,
                               byrow = FALSE))
    names(samps) <- c("T0", "Tm", "c", "sigma")
  }
  else{
    samps <- data.frame(matrix(c(T0, Tm, cc), ncol = 3, byrow=FALSE))
    names(samps) <- c("T0", "Tm", "c")
  }
  
  return(samps)
}

make.quad.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000), sig=TRUE){
  T0<-Tm<-qd<-sigma<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    T0<-c(T0, coda.samps[[i]][l1:l2,1])
    Tm<-c(Tm, coda.samps[[i]][l1:l2,2])
    qd<-c(qd, coda.samps[[i]][l1:l2,3])
    if(sig) sigma<-c(sigma, coda.samps[[i]][l1:l2,4])
  }
  if(sig){
    samps<-data.frame(matrix(c(T0, Tm, qd, sigma), ncol=4, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "qd", "sigma")
  }
  else{
    samps<-data.frame(matrix(c(T0, Tm, qd), ncol=3, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "qd")
    
  }
  return(samps)
}

briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
  }
  b
}
quad.2<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<--qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-0}
  }
  b
}

make.sims.temp.resp<-function(sim, samps, Temps, thinned, p.name="PDR", trunc.num=0){
  
  out<-data.sim<-NULL
  out<-list()
  data.sim<-matrix(NA, nrow=length(Temps), ncol=length(thinned))
  for(i in 1:length(thinned)){
    
    if(sim == "briere"){
      c<-as.numeric(samps$c[thinned[i]])
      Tm<-as.numeric(samps$Tm[thinned[i]])
      T0<-as.numeric(samps$T0[thinned[i]])
      w0<-which(Temps<=T0)
      wM<-which(Temps>=Tm)
      data.sim[,i]<-briere(Temps, c, Tm, T0)
      data.sim[c(w0,wM),i]<-0
    }
    
    if(sim == "briere.trunc"){
      c<-as.numeric(samps$c[thinned[i]])
      Tm<-as.numeric(samps$Tm[thinned[i]])
      T0<-as.numeric(samps$T0[thinned[i]])
      w0<-which(Temps<=T0)
      wM<-which(Temps>=Tm)
      data.sim[,i]<-briere.trunc(Temps, c, Tm, T0)
      data.sim[c(w0,wM),i]<-0
    }
    
    if(sim == "quad"){
      T0<-as.numeric(samps$T0[i])
      Tm<-as.numeric(samps$Tm[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad.2(Temps, T0=T0, Tm=Tm, qd=qd)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }
    
    if(sim == "quad.pos.trunc"){    # added by EAM on 8/31/15
      inter<-as.numeric(samps$inter[i])
      n.slope<- as.numeric(samps$n.slope[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad(Temps, inter=inter, n.slope=n.slope, qd=qd)
      w<-which(data.sim[,i]<trunc.num)
      data.sim[w,i]<-trunc.num
    }
    
    if(sim == "quad.trunc"){
      T0<-as.numeric(samps$T0[i])
      Tm<-as.numeric(samps$Tm[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad.trunc(Temps, T0=T0, Tm=Tm, qd=qd)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }
    
    if(sim == "linear"){
      n.inter<-as.numeric(samps$n.inter[i])
      slope<-as.numeric(samps$slope[i])
      data.sim[,i]<-linear(Temps, inter=-n.inter, slope=slope)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }
    
    if(sim == "nlinear"){
      inter<-as.numeric(samps$inter[i])
      n.slope<-as.numeric(samps$n.slope[i])
      data.sim[,i]<-linear(Temps, inter=inter, slope=-n.slope)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }
    
  }
  
  list(param = p.name,
       T   = Temps,
       fits = data.sim)
  
}

temp.sim.quants<-function(sim.data, l.temps, byCol=FALSE,
                          probs=c(0.025, 0.975)){
  
  q<-matrix(NA, nrow=length(probs), ncol=l.temps)
  if(byCol) for(i in 1:l.temps) q[,i]<-quantile(sim.data[,i], probs, na.rm=TRUE)
  else for(i in 1:l.temps) q[,i]<-quantile(sim.data[i,], probs, na.rm=TRUE)
  
  return(q)
}

# the following code is given by Villena et al. for fitting JAGS models of the
# diffferent types:

# Briere model (MDR, PDR, a)
jags_briere.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)T(0,)
mu.temp[i] <- c*T[i]*(T[i]-T0)*sqrt((Tm-T[i])*(Tm>T[i]))
mu[i] <- 0*(mu.temp[i]<0) + mu.temp[i]*(mu.temp[i]>0)

}

c ~ dgamma(1,10)
Tm ~ dunif(25,45)
T0  ~ dunif(0, 24)
sigma<-1/tau
tau ~ dgamma(0.0001, 0.0001)

}"


# Concave down quadratic model (PEA, EFD, bc)
jags_quad.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)
# Y[i] ~ dnorm(mu[i], tau)T(0,)
mu[i] <- -qd*(T[i]-T0)*(T[i]-Tm)*((T[i]>T0))*((T[i]<Tm))
}

Tm  ~ dunif(25,45)
T0 ~ dunif(0,24)
qd  ~ dgamma(1,1)
sigma<-1/tau
tau ~ dgamma(0.0001, 0.0001)

}"

define_jags_model <- function(data,
                              model_text,
                              mcmc_params,
                              inits) {
  jags.model(textConnection(model_text),
             data = list(
               Y = data$trait,
               T = data$T,
               N = length(data$T)
             ),
             n.chains = mcmc_params$n_chains,
             inits = inits,
             n.adapt = mcmc_params$n_adapt) 
  
}

# spline through the prediction and temperaures, and then constrain to be
# non-negative (as spline can induce negatives)
positive_spline <- function(pred, temps_out) {
  
  function_raw <- splinefun(temps_out, pred)
  function_positive <- function(temperature) {
    pmax(0, function_raw(temperature))
  }
  
  function_positive
  
}

# make a positive-constrained function from MCMC samples, for the required model
# type
make_function_jags <- function(coda_samples,
                          mcmc_params,
                          model_type = c("briere", "quad"),
                          temps_out = seq(-20, 80, by = 0.1)) {
  
  model_type <- match.arg(model_type)
  
  sample_function <- switch(model_type,
                            briere = make.briere.samps,
                            quad = make.quad.samps)
  
  # This command combines the samples from the chains into a format that we can
  # use for further analyses. Use appropriate model for specific traits
  samps <- sample_function(coda_samples,
                           nchains = mcmc_params$n_chains,
                           samp.lims = c(1, mcmc_params$n_samps))
  
  # use the parameter samples to get posterior samples of the temperature
  # response curve, and compute the posterior mean curve
  out <- make.sims.temp.resp(sim = model_type,
                             samps,
                             temps_out,
                             thinned = seq(1,
                                           mcmc_params$n_samps,
                                           length = 1000))
  
  post_mean <- rowMeans(out$fits)
  
  positive_spline(post_mean, temps_out)
  
}

# Use the Villena et al JAGS code to estimate the relationship between
# temperature and aquatic development rate for the required species
fit_mdr_temp <- function(data,
                         species = c("An. stephensi", "An. gambiae"),
                         # specify the parameters that control the MCMC
                         mcmc_params = list(
                           n_chains = 5,
                           n_adapt = 10000,
                           n_samps = 20000,
                           n_burn = 10000
                         ),
                         plot_fit = TRUE
) {
  species <- match.arg(species)
  
  data_sub <- data %>%
    filter(
      trait.name == "mdr",
      specie == species
    )
  
  # Use the corresponding JAGS model for each trait
  model <- define_jags_model(data_sub,
                             jags_briere.bug,
                             mcmc_params,
                             inits = list(
                               Tm = 31,
                               T0 = 5,
                               c = 0.00007
                             ))
  
  update(model, mcmc_params$n_samps)
  model_samps_coda <- coda.samples(model,
                                   c("c", "Tm", "T0", "sigma"),
                                   mcmc_params$n_samps)
  if (plot_fit) {
    # visual check for convergence
    plot(model_samps_coda, ask = TRUE) 
  }
  
  make_function_jags(model_samps_coda,
                     mcmc_params,
                     model_type = "briere")
  
}

# Use the Villena et al JAGS code to estimate the relationship between
# temperature and eggs per female per day for the required species

# note that we modell the EFD curve as a down quadratic, as described in the
# paper and SI, but the one plotted in the MS is clearly a Gaussian. I also had
# to remove the observation truncation in the model definition, since a bunch of
# 0s were observed and these were being thrown out, and also to flip the sign on
# the 'quad.2' function above to match the downward quadratic
fit_efd_temp <- function(data,
                         species = c("An. stephensi", "An. gambiae"),
                         # specify the parameters that control the MCMC
                         mcmc_params = list(
                           n_chains = 5,
                           n_adapt = 10000,
                           n_samps = 20000,
                           n_burn = 10000
                         ),
                         plot_fit = TRUE
) {
  species <- match.arg(species)
  
  data_sub <- data %>%
    filter(
      trait.name == "efd",
      specie == species
    )
  
  # do EFD as convex down
  model <- define_jags_model(data_sub,
                             jags_quad.bug,
                             mcmc_params,
                             inits = list(
                               Tm = 31,
                               T0 = 5,
                               qd = 0.00007
                             ))
  
  update(model, mcmc_params$n_samps)
  model_samps_coda <- coda.samples(model,
                                   c("qd", "Tm", "T0", "sigma"),
                                   mcmc_params$n_samps)
  if (plot_fit) {
    # visual check for convergence
    plot(model_samps_coda, ask = TRUE) 
  }
  
  make_function_jags(model_samps_coda,
                     mcmc_params,
                     model_type = "quad")
  
}

# given an mgcv mortality model, return the function for probability of survival
make_function_mgcv <- function(mortality_model,
                               temps_out = seq(-20, 80, by = 0.1)) {
  pred_df <- data.frame(temperature = temps_out)
  daily_preds <- 1 - predict(mortality_model, pred_df, type = "response")
  positive_spline(daily_preds, temps_out)
}

# DAS: Daily probability of survival during aquatic stages, computed from PEA:
# probability of surviving from egg to adult as a function of temperature,
# refitting the PEA data as a daily rate, using a cox proportional hazards
# model, with the exposure time given by the expected time to emergence from the
# MDR model. Ie. the daily survival probability may be lower at temperatures
# where the MDR is high, since they need to survive for less long. This enables
# us to model larval survival on a much shorter timestep
fit_das_temp <- function(data,
                         mdr_temp_fun,
                         species = c("An. stephensi", "An. gambiae"),
                         plot_fit = TRUE) {
  
  data_sub <- data %>%
    filter(
      trait.name == "e2a",
      specie == species,
      !is.na(initial)
    )
  
  data_survival <- data_sub %>%
    mutate(
      survived = round(initial * trait),
      died = initial - survived,
      exposure_period = 1 / mdr_temp_fun(temperature = T),
      exposure_offset = log(exposure_period) 
    ) %>%
    rename(
      temperature = T
    )
  
  # model daily survival probabilities at these temperatures via a proportional
  # hazards model, with the exposure period given by the MDR model
  mortality_model <- mgcv::gam(
    cbind(died, survived) ~ s(temperature),
    family = stats::binomial("cloglog"),
    offset = exposure_offset,
    data = data_survival,
    # enforce extra smoothing to make this consistent with the others 
    gamma = 20
  )
  
  if (plot_fit) {
    plot(mortality_model)
  }
  
  # function to return function of temperature
  make_function_mgcv(mortality_model)
  
} 

# given the daily survival rate in aquatic stages (DAS) and the aquatic stage
# development rate (MDR), return the probability of surviving from an egg to an
# adult at all (PEA)
make_pea_temp <- function(das_temp, mdr_temp) {
  function(temperature) {
    das_temp(temperature) ^ (1 / mdr_temp(temperature))
  }
}


# load the required R packages
library(rjags)
library(MASS)
library(tidyverse)
library(mgcv)

# Note: I had to get jags running on my M2 mac, so I did the following before
# loading rjags:
# Install jags at terminal with: `brew install jags`
# then install rjags pointing to this jags (modify path to wherever jags is,
# found with `which jags` at terminal)
# devtools::install_url("http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
#                       args="--configure-args='--with-jags-include=/opt/homebrew/bin/jags/include/JAGS        
#                                               --with-jags-lib=/opt/homebrew/bin/jags/lib'")

# load data, provided in the supplemental information to Villena et al., and clean/augment it
data_villena <- read.csv("data/life_history_params/oswaldov-Malaria_Temperature-16c9d29/data/traits.csv",
                     header = TRUE,
                     row.names = 1) %>%
  # Note, I can't find a study with this name and year that does aquatic stage
  # survival, only adult survival, so I am assuming this is a mistake and removing
  # it (6 observations)
  filter(
    !(trait.name == "e2a" & ref == "Murdock et al. 2016") 
  ) %>%
  # Ditto this study, there is a 2004 study by these authors which does adult
  # survival (3 observations) https://doi.org/10.1079/ber2004316
  filter(
    !(trait.name == "e2a" & ref == "Kirby and Lindsay 2009") 
  ) %>%
  # add on the initial number of eggs in the aquatic survival experiements (from
  # going back to the literature)
  mutate(
    initial = case_when(
      # https://doi.org/10.1111/gcb.12240
      ref == "Paaijmans et al. 2013" & trait.name == "e2a" ~ 50,
      # https://doi.org/10.1079/BER2003259
      ref == "Bayoh and Lindsay 2003" & trait.name == "e2a" ~ 30, 
      # J VBDs, no doi, initial number not given:
      # https://www.mrcindia.org/journal/issues/464295.pdf
      ref == "Olayemi and Ande 2009" & trait.name == "e2a" ~ NA,
      .default = NA
    )
  )

# use Villena et al. code to fit functions for MDR and EFD, against temperature
# for An. stephensi, and refit PEA ina. way that enables daily survival
# probabilities to be computed

# MDR: mosquito development rate (time to move through aquatic stages from egg
# to adult) as a function of temperature
mdr_temp_As <- fit_mdr_temp(data_villena, species = "An. stephensi")
mdr_temp_Ag <- fit_mdr_temp(data_villena, species = "An. gambiae")

# EFD: eggs per female per day as a function of temperature
efd_temp_As <- fit_efd_temp(data_villena, species = "An. stephensi")
efd_temp_Ag <- fit_efd_temp(data_villena, species = "An. gambiae")

# DAS: daily probability of survival in the aquatic stages (eggs, larvae, pupae)
das_temp_As <- fit_das_temp(data_villena,
                            mdr_temp_As,
                            species = "An. stephensi")
das_temp_Ag <- fit_das_temp(data_villena,
                            mdr_temp_As,
                            species = "An. gambiae")

# PEA: overall probability of surviving from an egg to an adult
pea_temp_As <- make_pea_temp(das_temp_As, mdr_temp_As)
pea_temp_Ag <- make_pea_temp(das_temp_Ag, mdr_temp_Ag)


# do for both species, and do this in ggplot, also with overlays between species

curves <- expand_grid(
  temperature = seq(0, 60, by = 1),
  trait_name = c("MDR", "PEA", "EFD"),
  species = c("An. stephensi", "An. gambiae")
) %>%
  mutate(
    trait = case_when(
      trait_name == "MDR" & species == "An. stephensi" ~ mdr_temp_As(temperature),
      trait_name == "MDR" & species == "An. gambiae" ~ mdr_temp_Ag(temperature),
      trait_name == "PEA" & species == "An. stephensi" ~ pea_temp_As(temperature),
      trait_name == "PEA" & species == "An. gambiae" ~ pea_temp_Ag(temperature),
      trait_name == "EFD" & species == "An. stephensi" ~ efd_temp_As(temperature),
      trait_name == "EFD" & species == "An. gambiae" ~ efd_temp_Ag(temperature),
    )
  )

data_villena_plot <- data_villena %>%
  mutate(
    trait_name = toupper(trait.name),
    trait_name = case_when(
      trait_name == "E2A" ~ "PEA",
      .default = trait_name
    ),
    species = specie,
    temperature = T
  ) %>%
  filter(
    trait_name %in% curves$trait_name,
    species %in% curves$species
  )

ggplot(
  curves,
  aes(y = trait,
      x = temperature)
) +
  geom_line() +
  facet_grid(trait_name ~ species,
             scales = "free_y") +
  theme_minimal() +
  geom_point(
    data = data_villena_plot,
    alpha = 0.3
  )


# now add density dependence effect to daily aquatic survival (from egg to adult)
# too - From Evans Figure 1D&F
library(tidyverse)
stephensi_survival <- read_csv("data/life_history_params/evans/data/clean/CSVs/survival.csv",
                               show_col_types = FALSE) %>%
  filter(
    Species == "Stephensi",
    AeDens == 0
  ) %>%
  mutate(
    NumInitial = StDens / 2
  ) %>%
  select(
    sex = Sex,
    temperature = Temp,
    density = StDens,
    replicate = Replicate,
    initial = NumInitial,
    survived = NumSurvived
  ) %>%
  mutate(
    # add a random ID for the experimental trial (M and F in together, but
    # multiple replicates per trial)
    trial = paste(temperature, density, replicate, sep = "_"),
    trial_id = as.numeric(factor(trial)),
    # add in the logit survival probability (from egg to adult) at no density as
    # a covariate for the temperature effect
    logit_survival_zero_density = qlogis(das_temp_As(temperature)),
    # remove the intercept term, and just add a parameter for the males
    is_male = as.numeric(sex == "Male")
  )

# model this
library(lme4)
stephensi_dd_model <- glmer(
  cbind(survived, initial - survived) ~ -1 +
    # no intercept, but a scaling factor on the Villena survival prob to get the
    # 0-density curve in terms of the survival prob used here (cumulative over
    # an undisclosed period)
    logit_survival_zero_density +
    # dummy for sex differences
    is_male +
    # logit-linear effect of density
    density +
    # random effect for the trial id (unique per replicate/condition
    # interaction)
    (1|trial_id),
  family = binomial,
  data = stephensi_survival
)

summary(stephensi_dd_model)

# residuals look fine (not surprising, given random effects)
plot(stephensi_dd_model)

# female is the default, and we don't care about the scaling on the survival
# probability (the daily egg to adult survival prob from Villena is more
# useful), so just scale that by the density term on the logit scale
density_effect <- fixef(stephensi_dd_model)[["density"]]

# create final function
das_temp_dens_As <- function(temperature, density) {
  temp_effect <- das_temp_As(temperature)
  plogis(qlogis(temp_effect) + density_effect * density)
}

# plot the model
survival_summary <- stephensi_survival %>%
  filter(
    sex == "Female"
  ) %>%
  group_by(temperature, density) %>%
  summarise(
    across(c(initial, survived), sum),
    .groups = "drop"
  ) %>%
  mutate(
    prob = survived / initial
  ) %>%
  select(-initial, -survived)

temp <- seq(10, 40, length.out = 100)
dens <- seq(0, 128, length.out = 100)
all <- expand_grid(dens, temp)
resp <- apply(all, 1, function(x) das_temp_dens_As(temperature = x["temp"],
                                               density = x["dens"]))
dim(resp) <- c(length(temp), length(dens))
contour(x = temp,
        y = dens,
        z = resp,
        levels = c(0.01, seq(0.1, 0.9, length.out = 9)),
        xlab = "Temperature (celsius)",
        ylab = "Density (larvae per 250ml)",
        main = "Daily probability of survival (egg to emergence)\nas a function of temperature and larval density")

points(density ~ temperature,
       data = survival_summary,
       pch = 21,
       bg = grey(1 - prob),
       col = grey(0.4))


# save the data for all of these functions for later use in modelling
mdr_function_data <- data.frame(temperature = Temps,
                                value = mdr_function(Temps))

efd_function_data <- data.frame(temperature = Temps,
                                value = efd_function(Temps))

das_function_data <- data.frame(temperature = Temps,
                                raw_value = das_temp_As(Temps),
                                density_coef = density_effect)

# can't save the functions as the environments get lost, so save the data and
# reconstitute the functions later

saveRDS(mdr_function_data, file = "data/life_history_params/mdr_function_data.RDS")
saveRDS(das_function_data, file = "data/life_history_params/das_function_data.RDS")
saveRDS(efd_function_data, file = "data/life_history_params/efd_function_data.RDS")
