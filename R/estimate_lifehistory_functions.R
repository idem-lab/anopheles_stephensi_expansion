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

# density dependence effects on daily aquatic survival (from egg to
# adult) for An stephensi from Evans figure 1 panels D&F
# https://doi.org/10.1002/eap.2334
load_stephensi_survival_data <- function(){
  
  read_csv("data/life_history_params/evans/data/clean/CSVs/survival.csv",
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
      # remove the intercept term, and just add a parameter for the males
      is_male = as.numeric(sex == "Male")
    )
  
}

# load data on temperatire-dependence of life history traits, provided in the
# supplemental information to Villena et al., and clean andaugment it
load_villena_data <- function() {
  
  read.csv("data/life_history_params/oswaldov-Malaria_Temperature-16c9d29/data/traits.csv",
           header = TRUE,
           row.names = 1) %>%
    # I can't find a study with this name and year that does aquatic stage
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
  
}

# given a survival temperature function and a density dependence log hazard
# coeficient, return a survival function of temperature and density
make_surv_temp_dens_function <- function(surv_temp_function, dd_effect) {
  function(temperature, density) {
    # get daily *mortality probability* at zero/low density at this temperature
    daily_mortality_zero_density <- 1 - surv_temp_function(temperature)
    # convert to the log hazard for a single day
    loghaz_mortality_zero_density <- log(-log(1 - daily_mortality_zero_density))
    # add on the density effect
    loghaz_mortality <- loghaz_mortality_zero_density + dd_effect * density
    # convert back to a daily *mortality probability*, including the density effect
    daily_mortality <- 1 - exp(-exp(loghaz_mortality))
    # and return as daily survival
    1 - daily_mortality
  }
}

dehydrate_lifehistory_function <- function(fun, path_to_object) {
  
  arguments <- formals(fun)
  e <- environment(fun)
  
  # determine how to store the required components
  if (identical(names(arguments), "temperature")) {
    object <- list(
      arguments = arguments,
      temperature_function_raw = e$function_raw,
      rectifier = function(x) {pmax(0, x)},
      dummy_function = function() {
        object$rectifier(
          object$temperature_function_raw(
            temperature
          )
        )
      }
    )
  } else if (identical(names(arguments), c("temperature", "density"))) {
    object <- list(
      arguments = arguments,
      temperature_function_raw = e$surv_temp_function,
      dd_effect = e$dd_effect,
      dummy_function = function() {
        daily_mortality_zero_density <- 1 - object$temperature_function_raw(temperature)
        loghaz_mortality_zero_density <- log(-log(1 - daily_mortality_zero_density))
        loghaz_mortality <- loghaz_mortality_zero_density + object$dd_effect * density
        daily_mortality <- 1 - exp(-exp(loghaz_mortality))
        1 - daily_mortality
      }
    )
  } else {
    stop("cannot dehydrate this function")
  }
  
  saveRDS(object, path_to_object)
  
}

# load An. gambiae adult survival data under temperature and humidity treatments
# from Bayoh's thesis
load_bayoh_data <- function() {
  read_csv(
    "data/life_history_params/adult_survival/bayoh/bayoh_an_gambiae_adult_survival.csv",
    col_types = cols(
      temperature = col_double(),
      humidity = col_double(),
      sex = col_character(),
      time = col_double(),
      died_cumulative = col_double(),
      alive = col_double()
    )) %>%
    mutate(
      sex = case_when(
        is.na(sex) ~ "mixed",
        .default = sex),
      replicate = 1,
      species = "An. gambiae",
      study = "bayoh"
    )
}

# load data from Krajacich et al. 2020
# https://doi.org/10.1186/s13071-020-04276-y on An gambiae adult survival under
# temperature and humidity combinations, when attempting to induce aestivation,
# and when not.
load_krajacich_data <- function() {
  krajacich <- read_csv("data/life_history_params/adult_survival/krajacich/aestivation.manu.files.scripts/22-Jan-18-R.formatted.masterlist.csv",
                        col_types = cols(
                          Experiment = col_character(),
                          Primed.as = col_character(),
                          Primed = col_character(),
                          Temp = col_character(),
                          `Date of death` = col_double(),
                          Censor = col_double()
                        )) %>%
    # temperatures and humidities were variable for some of these, so remove and
    # keep only the constant and clearly recorded temperatures
    filter(
      Temp != "SE",
      Temp != "18.male"
    ) %>%
    mutate(
      Temp = as.numeric(Temp)
    ) %>%
    # combine priming and experiments to get different replicates
    mutate(
      replicate = paste(Experiment, Primed.as, Primed, sep = "_"),
      replicate = match(replicate, unique(replicate))
    ) %>%
    select(
      -Experiment,
      -Primed.as,
      -Primed
    ) %>%
    rename(
      temperature = Temp,
      time = `Date of death`,
      status = Censor
    )
    

  # each row is an individual mosquito, Date of death is actually the last day
  # in the timeseries for that mosquito, and status is whether they were alive
  # (ie. when the experiment stopped) or dead (ie. day is the day they died)
  # then. We need to get the number alive at the start of each day (per
  # experiment and temp and humidity), and the number dying on that day
  
  # get all possible days, for all experiments
  all_combos <- expand_grid(
    replicate = unique(krajacich$replicate),
    temperature = unique(krajacich$temperature),
    time = seq(min(krajacich$time), max(krajacich$time))
  )
  
  # get the number of mosquitos at the start of each of these experiments
  starting <- krajacich %>%
    group_by(
      replicate,
      temperature
    ) %>%
    summarise(
      starting = n(),
      .groups = "drop"
    )
  
  # and the numbers dying on each day that one or more died on
  died <- krajacich %>%
    group_by(
      replicate,
      temperature,
      time
    ) %>%
    summarise(
      died = sum(status == 1),
      .groups = "drop"
    ) 
  
  # pull these all together, and add on other info
  all_combos %>%
    left_join(
      starting,
      by = join_by(replicate, temperature)
    ) %>%
    left_join(
      died,
      by = join_by(replicate, temperature, time)
    ) %>%
    mutate(
      died = replace_na(died, 0)
    ) %>%
    arrange(
      replicate,
      temperature,
      time
    ) %>%
    group_by(
      replicate,
      temperature
    ) %>%
    mutate(
      died_cumulative = cumsum(died)
    ) %>%
    ungroup() %>%
    mutate(
      alive = starting - died_cumulative
    ) %>%
    select(
      -starting,
      -died
    ) %>%
    mutate(
      sex = "F",
      species = "An. gambiae",
      humidity = 85,
      study = "krajacich"
    )
}

# load An stephensi adult survival data under temperature treatments from
# Miazgowicz et al. 2020 https://doi.org/10.1098/rspb.2020.1093, from Data
# dryad: https://doi.org/10.5061/dryad.8cz8w9gmd and prepare for modelling
load_miazgowicz_data <- function() {
  miazgowicz <- read_csv("data/life_history_params/adult_survival/miazgowicz/constant_master.csv",
                         col_types = cols(
                           Date = col_character(),
                           Time = col_time(format = ""),
                           Treatment = col_double(),
                           Block = col_double(),
                           Donor = col_double(),
                           Female = col_double(),
                           Feed = col_double(),
                           Size = col_character(),
                           Laid = col_double(),
                           Count = col_double(),
                           Dead = col_double(),
                           Day = col_double()
                         )) %>%
    group_by(
      Treatment,
      Block,
      Day
    ) %>%
    # Dead column is coded as 1 for dead, 0 for alive, and 2 for censored
    # (alive, but study discontinued), so recode these censored observations as alive, and count the number of observations
    # mutate(
    #   dead = as.numeric(Dead != 1)
    # ) %>%
    summarise(
      alive = n_distinct(Female[Dead %in% c(0, 2)]),
      died = n_distinct(Female[Dead %in% c(1)]),
      .groups = "drop"
    ) %>%
    group_by(
      Treatment,
      Block
    ) %>%
    mutate(
      died_cumulative = cumsum(died),
      humidity = 80,
      sex = "F",
      species = "An. stephensi",
      study = "miazgowicz"
    ) %>%
    ungroup() %>%
    select(-died) %>%
    rename(
      replicate = Block,
      temperature = Treatment,
      time = Day
    )
  
}

# load temperature-dependent data from Shapiro et al. 2017
# https://doi.org/10.1371/journal.pbio.2003489 from the data uploaded to data
# dryad https://doi.org/10.5061/dryad.74839
load_shapiro_data <- function() {
  
  shapiro <- read_csv("data/life_history_params/adult_survival/shapiro/temp.surv.csv",
                      col_types = cols(
                        id = col_double(),
                        expt = col_double(),
                        temp = col_double(),
                        cup = col_double(),
                        day = col_double(),
                        status = col_double()
                      )) %>%
    # combine experimental blocks and cups into a single replicate effect
    mutate(
      replicate = paste(expt, cup, sep = "_"),
      replicate = match(replicate, unique(replicate))
    ) %>%
    rename(
      temperature = temp,
      time = day
    )
  
  # each row is an individual mosquito, day is the last day in the timeseries
  # for that mosquito, and status is whether they were alive (ie. when the
  # experiment stopped) or dead (ie. day is the day they died) then. We need to
  # get the number alive at the start of each day (per experiment/cup and temp),
  # and the number dying on that day
  
  # get all possible days, for all experiments
  all_combos <- expand_grid(
    replicate = unique(shapiro$replicate),
    temperature = unique(shapiro$temperature),
    time = seq(min(shapiro$time), max(shapiro$time))
  )
  
  # get the number of mosquitos at the start of each of these experiments
  starting <- shapiro %>%
    group_by(
      replicate,
      temperature
    ) %>%
    summarise(
      starting = n(),
      .groups = "drop"
    )
  
  # and the numbers dying on each day that one or more died on
  died <- shapiro %>%
    group_by(
      replicate,
      temperature,
      time
    ) %>%
    summarise(
      died = sum(status == 1),
      .groups = "drop"
    ) 
  
  # pull these all together, and add on other info
  all_combos %>%
    left_join(
      starting,
      by = join_by(replicate, temperature)
    ) %>%
    left_join(
      died,
      by = join_by(replicate, temperature, time)
    ) %>%
    mutate(
      died = replace_na(died, 0)
    ) %>%
    arrange(
      replicate,
      temperature,
      time
    ) %>%
    group_by(
      replicate,
      temperature
    ) %>%
    mutate(
      died_cumulative = cumsum(died)
    ) %>%
    ungroup() %>%
    mutate(
      alive = starting - died_cumulative
    ) %>%
    select(
      -starting,
      -died
    ) %>%
    mutate(
      sex = "F",
      species = "An. stephensi",
      humidity = 80,
      study = "shapiro"
    )
    
}



# load the required R packages

# for modelling
library(rjags)
library(mgcv)
library(lme4)
library(MASS)

# for data handling and plotting
library(tidyverse)

# for contour labelling
library(metR)


# Note: I had to get jags running on my M2 mac, so I did the following before
# loading rjags:
# Install jags at terminal with: `brew install jags`
# then install rjags pointing to this jags (modify path to wherever jags is,
# found with `which jags` at terminal)
# devtools::install_url("http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
#                       args="--configure-args='--with-jags-include=/opt/homebrew/bin/jags/include/JAGS        
#                                               --with-jags-lib=/opt/homebrew/bin/jags/lib'")

# load data, provided in the supplemental information to Villena et al., and clean/augment it
data_villena <- load_villena_data()

# use Villena et al. code to fit functions for MDR and EFD, against temperature
# for An. stephensi, and refit PEA in a way that enables daily survival
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
                            mdr_temp_fun = mdr_temp_As,
                            species = "An. stephensi")
das_temp_Ag <- fit_das_temp(data_villena,
                            mdr_temp_fun = mdr_temp_As,
                            species = "An. gambiae")

# PEA: overall probability of surviving from an egg to an adult (used only for
# plotting)
pea_temp_As <- make_pea_temp(das_temp_As, mdr_temp_As)
pea_temp_Ag <- make_pea_temp(das_temp_Ag, mdr_temp_Ag)

# plot fitted curves and data
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
    ),
    trait_name = case_when(
      trait_name == "MDR" ~ "Aquatic development",
      trait_name == "PEA" ~ "Aquatic survival",
      trait_name == "EFD" ~ "Adult eggs/day"
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
    temperature = T,
    trait_name = case_when(
      trait_name == "MDR" ~ "Aquatic development",
      trait_name == "PEA" ~ "Aquatic survival",
      trait_name == "EFD" ~ "Adult eggs/day"
    )
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
             scales = "free_y",
             switch = "y") +
  geom_point(
    data = data_villena_plot,
    alpha = 0.3
  ) +
  theme_minimal() +
  theme(strip.placement = "outside") +
  xlab("Temperature (C)") +
  ylab("")

ggsave("figures/lifehistory_temperature.png",
       bg = "white",
       width = 5,
       height = 5)

ggplot(
  filter(curves,
         species == "An. stephensi"),
  aes(y = trait,
      x = temperature)
) +
  geom_line() +
  facet_grid(trait_name ~ species,
             scales = "free_y",
             switch = "y") +
  geom_point(
    data = filter(data_villena_plot,
                  species == "An. stephensi"),
    alpha = 0.3
  ) +
  theme_minimal() +
  theme(strip.placement = "outside") +
  xlab("Temperature (C)") +
  ylab("")

ggsave("figures/lifehistory_temperature_stephensi.png",
       bg = "white",
       width = 3,
       height = 5)

ggplot(
  filter(curves,
         species == "An. gambiae"),
  aes(y = trait,
      x = temperature)
) +
  geom_line() +
  facet_grid(trait_name ~ species,
             scales = "free_y",
             switch = "y") +
  geom_point(
    data = filter(data_villena_plot,
                  species == "An. gambiae"),
    alpha = 0.3
  ) +
  theme_minimal() +
  theme(strip.placement = "outside") +
  xlab("Temperature (C)") +
  ylab("")

ggsave("figures/lifehistory_temperature_gambiae.png",
       bg = "white",
       width = 3,
       height = 5)

# now add density dependence effect to daily aquatic survival (from egg to
# adult) too for An stephensi

# load density and temperature treatment survival data from Evans figure 1
# panels D&F and use it to compute a density coefficent on the log hazard scale
# for the probability of survival. Note that Evans don't give the data in
# survival curve form and don't specify the durations of the experiments so we
# again use the modelled MDR as the expoure period and a Cox PH model to model
# that. Note that this assumes density does not affect development rate, but any
# impact of that on fraction reaching adulthood will be accounted for in the
# survival effect. Evan et al.: https://doi.org/10.1002/eap.2334
stephensi_survival <- load_stephensi_survival_data() %>%
  mutate(
    died = initial - survived,
    exposure_period = 1 / mdr_temp_As(temperature = temperature),
    exposure_offset = log(exposure_period),
    mortality_zero_density = 1 - das_temp_As(temperature),
    cloglog_mortality_zero_density = log(-log(1 - mortality_zero_density))
  )

mortality_model <- lme4::glmer(
  cbind(died, survived) ~ 1 +
    # dummy for sex differences
    is_male +
    # logit-linear effect of density
    density +
    # random effect for the trial id (unique per replicate/condition
    # interaction)
    (1|trial_id),
  family = stats::binomial("cloglog"),
  offset = exposure_offset + cloglog_mortality_zero_density,
  data = stephensi_survival
)

# residuals look fine (not surprising, given all the random effects)
plot(mortality_model)

# this is the density effect on *mortality* in the COX PH model
# female is the default, and we don't care about the scaling on the survival
# probability (the daily egg to adult survival prob from Villena is more
# useful), so just scale that by the density term on the logit scale
dd_effect_As <- fixef(mortality_model)[["density"]]

# Note that we could roll all of this into the GAM temperature survival model,
# but that would be hard to model with multiple datasts and noise terms, and
# would require a 2D spline which would be more computationally costly to
# predict from when simulating the model

# create the final function
das_temp_dens_As <- make_surv_temp_dens_function(
  surv_temp_function = das_temp_As,
  dd_effect = dd_effect_As
)

# plot the survival model
density_plotting <- expand_grid(
  temperature = seq(10, 40, length.out = 100),
  density = seq(0, 210, length.out = 100)  
) %>%
  rowwise() %>%
  mutate(
    prob = das_temp_dens_As(temperature, density)
  )

max_grey <- 0.6
cols <- grey(1 - max_grey * seq(0, 1, by = 0.1))
cols[1] <- "transparent"

ggplot(density_plotting,
       aes(y = density,
           x = temperature,
           z = prob)) +
  geom_contour_filled(
    binwidth = 0.1
  ) +
  geom_contour(
    binwidth = 0.1,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    binwidth = 0.1,
    nudge_y = -5,
    # rotate = FALSE,
    skip = 0
    # label.placer = label_placer_n(3)
  ) +
  scale_fill_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Daily survival probability of aquatic stages") +
  ylab("Individuals per 250ml") +
  xlab("Temperature (C)")

ggsave("figures/aquatic_survival_stephensi.png",
       bg = "white",
       width = 5,
       height = 5)

# model adult survival as a function of air temperature and humidity. Reanalyse
# Bayoh data on an gambiae, include other data on longer-lived An gambiae, and
# jointly model An stephensi (with temperature effects, but fixed humidity.)


# fit temperature- and humidity-dependent survival curve

# load Mohammed Bayoh's thesis data for An gambiae with temperature and humidity
# treatments, 16CSS strain originally colonised from Lagos, Nigeria, in 1974
bayoh_Ag <- load_bayoh_data()

# load Krajacich et al 2020's data on inducing aestivation in An gambiae (a
# younger colony than Bayoh; coluzzi from Mali in 2012 and ss from Cameroon in
# 2008) strangely, the species is not distinguished in the results or data
# https://doi.org/10.1186/s13071-020-04276-y
krajacich_Ag <- load_krajacich_data()

# Villena has two datasets for temperature dependent An stephensi adult
# mortality rates, Shapiro et al. 2017
# https://doi.org/10.1371/journal.pbio.2003489 and Miazgowicz et al. 2020
# https://doi.org/10.1098/rspb.2020.1093 (listed as Kerri 2019, presumably
# referring to the preprint). The only adult survival data is from Bayoh's
# thesis, which we already have.
# from the 'long-standing colony' at Walt Reed.
shapiro_As <- load_shapiro_data()

# load data from Miazgowicz et al. 2020 https://doi.org/10.1098/rspb.2020.1093,
# from Data dryad: https://doi.org/10.5061/dryad.8cz8w9gmd from a 'long-standing
# colony (~40 years) of An. stephensi mosquitoes from Pennsylvania State
# University which were originally obtained from the Walter Reed Army Institute
# of Research'
miazgowicz_As <- load_miazgowicz_data()

# fit a proportional hazards model on probability of *mortality* via cloglog 
# trick: offset handles the cumulative effect, first term is baseline hazard, 
# temperature and humidity have a multiplicative effect, fixed intercepts for 
# each sex group

adult_survival_data <- bind_rows(
  bayoh_Ag,
  krajacich_Ag,
  miazgowicz_As,
  shapiro_As
) %>%
  mutate(
    species = factor(species),
    study = factor(study),
    # Define a preferred study for each species, to account for study
    # differences, without confounding the species differences. For An gambiae,
    # the preferred study is the one with the youngest colony. For An. stephensi
    # they are from the same colony (v. old), so just picking one at random
    non_preferred = case_when(
      study %in% c("krajacich", "miazgowicz") ~ 0,
      .default = 1
    ),
    id = row_number()
  ) %>%
  # treat each observation period as an independent observation, so account for
  # the number starting each time period (alive_start) and compute the number
  # that died (died_end) and number that lived (alive_end) in that period
  group_by(temperature, humidity, sex, replicate, species, study) %>%
  mutate(
    alive_start = c(0, alive[-n()]),
    died_end = c(0, diff(died_cumulative)),
    alive_end = alive_start - died_end,
    # how long was this survival interval?
    duration = c(0, diff(time)),
    .after = alive
  ) %>%
  ungroup() %>%
  # use the survival interval times the number of trials as the offset (NB approximation to the betabinomial)
  mutate(
    offset = log(duration * alive_start),
    As_mask = as.numeric(species == "An. stephensi"),
    temperature_As = temperature * As_mask
  ) %>%
  # remove rows for time zero starting values (anywhere there weren't mossies
  # alive at the start)
  filter(
    alive_start != 0,
    sex == "F"
  )

# we could model this as binomial with cloglog (survival model with proportional
# hazards), but there is additional dispersion even after fitting smooth terms,
# so we need extra observation-level variance. Betabinomial is not possible in.
# mgcv, so we use a negative binomial approximation to the BB

m <- mgcv::gam(died_end ~ 1 +
                 # intercept for species, and for the preferred study type (to
                 # account for differences in colony age)
                 species +
                 # non_preferred +
                 # nonlinear interaction between them
                 s(temperature, log(humidity), k = 15) +
                 # linear effect of mosquito age (duration of exposure)
                 time * study,
               # penalise all effects to zero (lasso-like regulariser)
               select = TRUE,
               data = adult_survival_data,
               # apply extra smoothing, to ad-hoc account for additional noise
               # in the observations
               gamma = 10,
               # fit with REML
               method = "REML",
               # account for the differing durations of the observation periods
               # using a cloglog and offset (proportional hazards model with
               # double censoring)
               offset = adult_survival_data$offset,
               family = mgcv::nb(link = "log")
               # set the scale estimator to avoid a bug (there shouldn't be a
               # scale parameter anyway)
               # control = list(scale.est = "deviance")
)
summary(m)
gratia::draw(m)
mgcv::gam.check(m)

# use randomised quantile residuals to check model fit

library(DHARMa)
sims <- simulateResiduals(m)
testOutliers(sims, type = "bootstrap")
plot(sims)

# mostly it's fine, but there are still some outliers

# get target-normally-distributed residuals, to explore lack of fit
resid_checks <- adult_survival_data %>%
  select(
    time,
    temperature,
    humidity,
    species,
    study
  ) %>%
  mutate(
    resid = qnorm(sims$scaledResiduals)
  )

# this should look uniform
hist(pnorm(resid_checks$resid))

# what proportion of datapoints are outliers?
round(100 * mean(!is.finite(resid_checks$resid)), 1)

# plot where these fall in time vs temp
ggplot(resid_checks,
       aes(y = jitter(time),
           x = jitter(temperature),
           color = is.finite(resid),
           size = !is.finite(resid),
           group = study)) +
  geom_point(
    alpha = 0.3,
    pch = 16
  ) +
  facet_grid(~ study)


# This is in a lab; field mortality is much greater. Charlwood (1997) estimates
# daily adult survival for An gambiae at around 0.83 in Ifakara with temperature
# = 25.6, humidity = 80%. Compute a correction on the hazard scale for the
# lab-based survival estimates.
mortality_prob_to_log_hazard <- function(mortality_prob) log(-log(1 - mortality_prob))
log_hazard_to_mortality_prob <- function(log_hazard) 1 - exp(-exp(log_hazard))

df_ifakara <- data.frame(temperature = 25.6,
                         humidity = 80,
                         sex = "F",
                         time = sqrt(.Machine$double.eps),
                         species = "An. gambiae",
                         temperature_As = 0,
                         replicate = 1,
                         id = 1,
                         non_preferred = 0,
                         # freshest An. gambiae population we have data for
                         study = "krajacich",
                         off = 0)

lab_daily_log_hazard <- predict(m, df_ifakara, type = "link")[1]
field_daily_log_hazard <- mortality_prob_to_log_hazard(1 - 0.83)
log_hazard_correction <- field_daily_log_hazard - lab_daily_log_hazard

# construct function of temperature and humidity
ds_temp_humid <- function(temperature, humidity, species, epsilon = sqrt(.Machine$double.eps)) {
  df <- data.frame(
    temperature = pmax(epsilon, temperature),
    humidity = pmax(epsilon, humidity),
    time = epsilon,
    sex = "F",
    id = 1,
    species = species,
    non_preferred = 0,
    study = ifelse(species == "An. gambiae",
                   "krajacich",
                   "miazgowicz"),
    temperature_As = case_when(
      species == "An. stephensi" ~ temperature,
      .default = 0
    ),
    replicate = 1
  )
  # browser()
  # lab_mortality_link <- predict(m, df, type = "link")
  # lab_mortality_prob <- 1 - exp(-exp(lab_mortality_link))
  lab_mortality_prob <- predict(m, df, type = "response")
  lab_daily_log_hazard <- mortality_prob_to_log_hazard(lab_mortality_prob)
  field_daily_log_hazard <- lab_daily_log_hazard + log_hazard_correction
  field_mortality_prob <- log_hazard_to_mortality_prob(field_daily_log_hazard)
  1 - field_mortality_prob
}

ds_temp_humid_Ag <- function(temperature, humidity) {
  ds_temp_humid(temperature, humidity, species = "An. gambiae")
}

ds_temp_humid_As <- function(temperature, humidity) {
  ds_temp_humid(temperature, humidity, species = "An. stephensi")
}

# plot the survival model

# contour bin width and colours
bw <- 0.05
max_grey <- 0.6
cols <- grey(1 - max_grey * seq(0, 1, by = bw))
cols[1] <- "transparent"

density_plotting <- expand_grid(
  temperature = seq(5, 40, length.out = 100),
  humidity = seq(0, 100, length.out = 100),
  species = c("An. stephensi", "An. gambiae")
) %>%
  mutate(
    prob = case_when(
      species == "An. stephensi" ~ ds_temp_humid_As(temperature, humidity),
      species == "An. gambiae" ~ ds_temp_humid_Ag(temperature, humidity),
      .default = NA
    ),
    prob_7d = prob ^ 7
  )


# summarise survivial data to overplot
survival_summary <- adult_survival_data %>%
  group_by(
    species,
    temperature,
    humidity,
    replicate,
    sex
  ) %>%
  ungroup() %>%
  group_by(
    species,
    temperature,
    humidity
  ) %>%
  summarise(
    alive = sum(alive_end),
    total = sum(alive_start),
    .groups = "drop"
  ) %>%
  mutate(
    survival = alive / total,
    survival_7d = survival ^ 7,
    survival_7d_cut = cut(survival_7d, seq(0, 1, by = bw))
  )

  
ggplot(density_plotting,
       aes(y = humidity,
           x = temperature)) +
  facet_wrap(~species) +
  geom_contour_filled(
    aes(z = prob_7d),
    binwidth = bw
  ) +
  geom_contour(
    aes(z = prob_7d),
    binwidth = bw,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    aes(z = prob_7d),
    binwidth = bw,
    nudge_y = -5,
    skip = 0
  ) +
  geom_point(
    data = survival_summary,
    mapping = aes(
      colour = survival_7d_cut
    ),
    shape = 16,
    size = 4
  ) +
  geom_point(
    data = survival_summary,
    shape = 21,
    colour = grey(0.4),
    size = 4
  ) +
  scale_fill_discrete(type = cols) +
  scale_color_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Adult survival vs air temperature and humidity",
          "Probability of surviving 1 week") +
  ylab("Humidity(%)") +
  xlab("Temperature (C)")

ggsave("figures/lifehistory_adult_survival.png",
       bg = "white",
       width = 7,
       height = 5)
# and export the survival model as before




# now we need to save these neatly to be loaded reused later, getting around all
# the awkward lexical scoping stuff. Use dehydrate_lifehistory_function() and
# rehydrate_lifehistory_function() to store and rejuvenate the functions as RDS
# files, with all relevant objects included

storage_path <- "data/life_history_params/dehydrated"

# MDR, EFD, DAS, PEA ~ temp
dehydrate_lifehistory_function(mdr_temp_As, file.path(storage_path, "mdr_temp_As.RDS"))
dehydrate_lifehistory_function(mdr_temp_Ag, file.path(storage_path, "mdr_temp_Ag.RDS"))
dehydrate_lifehistory_function(efd_temp_As, file.path(storage_path, "efd_temp_As.RDS"))
dehydrate_lifehistory_function(efd_temp_Ag, file.path(storage_path, "efd_temp_Ag.RDS"))
dehydrate_lifehistory_function(das_temp_As, file.path(storage_path, "das_temp_As.RDS"))
dehydrate_lifehistory_function(das_temp_Ag, file.path(storage_path, "das_temp_Ag.RDS"))
dehydrate_lifehistory_function(pea_temp_As, file.path(storage_path, "pea_temp_As.RDS"))
dehydrate_lifehistory_function(pea_temp_Ag, file.path(storage_path, "pea_temp_Ag.RDS"))

# DAS ~ temp + density
dehydrate_lifehistory_function(das_temp_dens_As, file.path(storage_path, "das_temp_dens_As.RDS"))

# use the following function to rehydrate them, e.g.:
#   storage_path <- "data/life_history_params/dehydrated"
#   das_temp_dens_As <- rehydrate_lifehistory_function(file.path(storage_path, "das_temp_dens_As.RDS"))

rehydrate_lifehistory_function <- function(path_to_object) {
  object <- readRDS(path_to_object)
  do.call(`function`,
          list(object$arguments,
               body(object$dummy_function)))
}

