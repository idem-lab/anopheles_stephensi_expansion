# Refit the models in villeno et al to get and interpolate the posterior
# predicted relationship against temperature of key life history parameters.
# Using the aadapted JAGS code from Villeno et al.

# I had to get jags running on my M1 mac, so I did the following:
# install jags at terminal with: `brew install jags` then install rjags pointing
# to this jags (modify path to wherever jags is, found with `which jags` at
# terminal)
# devtools::install_url("http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
#                       args="--configure-args='--with-jags-include=/opt/homebrew/bin/jags/include/JAGS        
#                                               --with-jags-lib=/opt/homebrew/bin/jags/lib'")

## First load the necessary R packages and files for the computation and
## visualization
library(rjags)
library(MASS)
library(tidyverse)

# grab this thing from here: https://rdrr.io/github/lorecatta/DENVclimate/src/R/mcmc_utils_all.R
# which seems not to exist on GH any more?
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

## Next load data
data.all <- read.csv("data/life_history_params/oswaldov-Malaria_Temperature-16c9d29/data/traits.csv",
                     header = TRUE,
                     row.names = 1)

## Select the trait of interest. Here we showed as an example MDR for An. gambiae
data_mdr <- data.all %>%
  filter(
    trait.name == "mdr",
    specie == "An. stephensi"
  )

data_pea <- data.all %>%
  filter(
    trait.name == "e2a",
    specie == "An. stephensi",
    # can't find a study with this name and year that does larval survival, only
    # adult
    ref != "Murdock et al. 2016" 
  )

data_efd <- data.all %>%
  filter(
    trait.name == "efd",
    specie == "An. stephensi"
  )

# # visualize your data
# plot(data_mdr$T, data_mdr$trait)
# plot(data_pea$T, data_pea$trait)

## specify the parameters that control the MCMC
n.chains <- 5
n.adapt <- 10000
n.samps <- 20000
n.burn <- 10000


##Choose the appropiate model for the specific trait 

##Briere model (MDR, PDR, a)

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


## Concave down quadratic model (PEA, EFD, bc)

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



## Use jags.model for the specific model with the appropiate
## default priors
mdr_model <- jags.model(textConnection(jags_briere.bug),
                    data = list(
                      Y = data_mdr$trait,
                      T = data_mdr$T,
                      N = length(data_mdr$T)
                    ),
                    n.chains = n.chains,
                    inits = list(
                      Tm = 31,
                      T0 = 5,
                      c = 0.00007
                    ),
                    n.adapt = n.adapt) 
update(mdr_model, n.samps)

mdr_model_samps_coda <- coda.samples(mdr_model,
                            c('c','Tm', 'T0', 'sigma'),
                            n.samps)
# ## check for convergence
# plot(mdr_model_samps_coda, ask = TRUE)

## This command combines the samples from the n.chains into a format
## that we can use for further analyses. Use appropiate model for specific traits
mdr_samps <- make.briere.samps(mdr_model_samps_coda,
                           nchains = n.chains,
                           samp.lims = c(1, n.samps))



# do EFD as convex down
efd_model <- jags.model(textConnection(jags_quad.bug),
                        data = list(
                          Y = data_efd$trait,
                          T = data_efd$T,
                          N = length(data_efd$T)
                        ),
                        n.chains = n.chains,
                        inits = list(
                          Tm = 31,
                          T0 = 5,
                          qd = 0.00007
                        ),
                        n.adapt = n.adapt) 
update(efd_model, n.samps)

efd_model_samps_coda <- coda.samples(efd_model,
                                     c('qd','Tm', 'T0', 'sigma'),
                                     n.samps)
# ## check for convergence
# plot(efd_model_samps_coda, ask = TRUE)

## This command combines the samples from the n.chains into a format
## that we can use for further analyses. Use appropiate model for specific traits
efd_samps <- make.quad.samps(efd_model_samps_coda,
                             nchains = n.chains,
                             samp.lims = c(1, n.samps))

## Next we want to use the parameter samples to get posterior samples
## of the temperature rsponses themselves
Temps <- seq(-20, 80, by = 0.1)
mdr_out <- make.sims.temp.resp(sim = "briere",
                               mdr_samps,
                               Temps,
                               thinned = seq(1, n.samps, length = 1000)) ## Example for MDR
efd_out <- make.sims.temp.resp(sim = "quad",
                               efd_samps,
                               Temps,
                               thinned = seq(1, n.samps, length = 1000)) ## Example for MDR

mdr_post_mean <- rowMeans(mdr_out$fits)
efd_post_mean <- rowMeans(efd_out$fits)

# spline these functions
mdr_function_raw <- splinefun(Temps, mdr_post_mean)
mdr_function <- function(temperature) {
  pmax(0, mdr_function_raw(temperature))
}

efd_function_raw <- splinefun(Temps, efd_post_mean)
efd_function <- function(temperature) {
  pmax(0, efd_function_raw(temperature))
}

plot(mdr_function, xlim = c(-20, 80), type = "l")
points(data_mdr$T, data_mdr$trait)
min(mdr_function(Temps))

plot(efd_function, xlim = c(-20, 80), type = "l", ylim = range(data_efd$trait))
points(data_efd$T, data_efd$trait)
min(efd_function(Temps))



# note that the fitted EFD curve is a down quadratic, as described in the paper
# and SI, but the one plotted in the MS is clearly a Gaussian. I also had to
# remove the observation truncation, since a bunch of 0s were observed and these
# were being thrown out, and flip the sign on the 'quad.2' function bove to
# match the downward quadratic


# refit the PEA data as a daily rate, using a cox proportional hazards model,
# with the exposure time given by the expected time to emergence from the MDR
# model. Ie. the daily survival probability can be lower at temperatures where
# the MDR is high, since they need to survive for less long. This enables us to
# model larval survival on a much shorter timestep

data_pea_survival <- data_pea %>%
  mutate(
    initial = 50,  # from Paaijmans
    survived = round(initial * trait),
    died = initial - survived,
    exposure_period = 1 / mdr_function(temperature = T),
    exposure_offset = log(exposure_period) 
  ) %>%
  rename(
    temperature = T
  )

library(mgcv)
# model daily survival probabilities at these temperatures via a proportional
# hazards model, with the exposure period given by the MDR model
pea_mortality_model <- gam(
  cbind(died, survived) ~ s(temperature),
  family = stats::binomial("cloglog"),
  offset = exposure_offset,
  data = data_pea_survival,
  # enforce extra smoothing to make this consistent with the others 
  gamma = 30
)

# plot(pea_mortality_model)

# predict daily survival to the temperature range
pred_df <- data.frame(temperature = Temps)
pea_daily_preds <- 1 - predict(pea_mortality_model, pred_df, type = "response")

# plot(pea_daily_preds ~ Temps, type = "l")
# rug(data_pea$T)

# # check it matches observed full-period mortalities by recombining with MDR
# pea_overall_preds <- pea_daily_preds ^ (1 / mdr_function(pred_df$temperature))
# 
# plot(pea_overall_preds ~ Temps, type = "l",
#      xlim = c(10, 40),
#      ylim = c(0, 1))
# points(data_pea$trait ~ data_pea$T)
# # good fit!


pea_function_temperature_raw <- splinefun(Temps, pea_daily_preds)
pea_function_temperature <- function(temperature) {
  pmax(0, pea_function_temperature_raw(temperature))
}
plot(pea_function_temperature, xlim = c(-20, 80), type = "l")
min(pea_function_temperature(Temps))

# now add density dependence effect to daily pea (survival from egg to adult)
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
    logit_survival_zero_density = qlogis(pea_function_temperature(temperature)),
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
pea_function <- function(temperature, density) {
  temp_effect <- pea_function_temperature(temperature)
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
resp <- apply(all, 1, function(x) pea_function(temperature = x["temp"],
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

pea_function_data <- data.frame(temperature = Temps,
                                raw_value = pea_function_temperature(Temps),
                                density_coef = density_effect)

# can't save the functions as the environments get lost, so save the data and
# reconstitute the functions later

saveRDS(mdr_function_data, file = "data/life_history_params/mdr_function_data.RDS")
saveRDS(pea_function_data, file = "data/life_history_params/pea_function_data.RDS")
saveRDS(efd_function_data, file = "data/life_history_params/efd_function_data.RDS")
