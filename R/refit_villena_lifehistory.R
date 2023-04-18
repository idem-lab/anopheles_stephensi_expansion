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
    data.all$trait.name == "mdr",
    data.all$specie == "An. stephensi"
  )

data_pea <- data.all %>%
  filter(
    data.all$trait.name == "e2a",
    data.all$specie == "An. stephensi"
  )

data_efd <- data.all %>%
  filter(
    data.all$trait.name == "efd",
    data.all$specie == "An. stephensi"
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
Y[i] ~ dnorm(mu[i], tau)T(0,)
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



# do PEA and EFD as convex down

pea_model <- jags.model(textConnection(jags_quad.bug),
                        data = list(
                          Y = data_pea$trait,
                          T = data_pea$T,
                          N = length(data_pea$T)
                        ),
                        n.chains = n.chains,
                        inits = list(
                          Tm = 31,
                          T0 = 5,
                          qd = 0.00007
                        ),
                        n.adapt = n.adapt) 
update(pea_model, n.samps)

pea_model_samps_coda <- coda.samples(pea_model,
                                     c('qd','Tm', 'T0', 'sigma'),
                                     n.samps)
# ## check for convergence
# plot(pea_model_samps_coda, ask = TRUE)

## This command combines the samples from the n.chains into a format
## that we can use for further analyses. Use appropiate model for specific traits
pea_samps <- make.quad.samps(pea_model_samps_coda,
                             nchains = n.chains,
                             samp.lims = c(1, n.samps))

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
update(efd_model, n.samps)

efd_model_samps_coda <- coda.samples(efd_model,
                                     c('qd','Tm', 'T0', 'sigma'),
                                     n.samps)
# ## check for convergence
plot(efd_model_samps_coda, ask = TRUE)

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
pea_out <- make.sims.temp.resp(sim = "quad",
                               pea_samps,
                               Temps,
                               thinned = seq(1, n.samps, length = 1000)) ## Example for MDR
efd_out <- make.sims.temp.resp(sim = "quad",
                               efd_samps,
                               Temps,
                               thinned = seq(1, n.samps, length = 1000)) ## Example for MDR

mdr_post_mean <- rowMeans(mdr_out$fits)
pea_post_mean <- rowMeans(pea_out$fits)
efd_post_mean <- rowMeans(efd_out$fits)

mdr_function_raw <- splinefun(Temps, mdr_post_mean)
mdr_function <- function(temperature) {
  pmax(0, mdr_function_raw(temperature))
}

pea_function_raw <- splinefun(Temps, pea_post_mean)
pea_function <- function(temperature) {
  pmax(0, pea_function_raw(temperature))
}

efd_function_raw <- splinefun(Temps, efd_post_mean)
efd_function <- function(temperature) {
  pmax(0, efd_function_raw(temperature))
}

plot(mdr_function, xlim = c(-20, 80), type = "l")
points(data_mdr$T, data_mdr$trait)
min(mdr_function(Temps))

plot(pea_function, xlim = c(-20, 80), type = "l")
points(data_pea$T, data_pea$trait)
min(pea_function(Temps))

plot(efd_function, xlim = c(-20, 80), type = "l", ylim = range(data_efd$trait))
points(data_efd$T, data_efd$trait)
min(efd_function(Temps))


saveRDS(mdr_function, file = "data/life_history_params/mdr_function.RDS")
saveRDS(pea_function, file = "data/life_history_params/pea_function.RDS")
saveRDS(efd_function, file = "data/life_history_params/efd_function.RDS")

