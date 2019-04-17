# Garbowski - 2017 CPNPP Traits
# Questions for Daniel Laughlin 
# March 30, 2019 

# Goals/What I'd appreciate help with (1 -3, maybe 4 down the road).
#Technical Questions:
  
# 1. How can I model variance (sigma or tau) from each population in JAGS? I FIGURED IT OUT! In "24-day model"
# 2. How can I include sampling date (harvest time, "H_num"" or "Days" in the data file) in my model. I'd like to compare estimates within/among population(s) but at different times? 
# 3. Is the model below appropriate for all traits? Are the priors/distributions reasonable? Even for traits with values 0-1?
# 4. Evenutally...How to compare among species?
  
# Below is an example of my process so far (translated to the best of my ability from STAN to JAGS).
# I am still "hacking" to try to get things to work. 

# Load packages
library(ggplot2)
library(rjags)
library(parallel)
library(runjags)
library(reshape)
library(coda)

# Load data 

Trait2017_data<-read.csv("/Users/MagdaGarbowski/Traits_2017_Grasses/Garbowski_CPNPP2017_Laughlin_Help/Garbowski_HECO_ELTR_PLPA_VUOC.csv")
# Subset Heterostipa comata (HECO) 
HECO_data<-Trait2017_data[Trait2017_data$SPECIES=="HECO",]
# Create integer index for POP_NUM
HECO_data$POP_ID.index<-as.integer(
factor(HECO_data[,c("POP_ID")],levels=c("HECO_AZNC_KV","HECO_UTC_KRR","HECO_UTC_TMR","HECO_UTEC_GR","HECO_UTSW_SP")))
HECO_data<-HECO_data[complete.cases(HECO_data), ]

#
#
#
#
#
#
#  Pooled Model - Pooled data for all populations of HECO. Working with SLA column.

# Input Data
data = list(
y.SLA = (HECO_data$SLA.TOT),
n = length(HECO_data$SLA.TOT)) # number of observations (n) for the JAGS counter

# Model
{
  sink("Garbowski_Pooled_POPS_JAGS.R")
  cat("
  model{
  # Priors
  mu_log_y ~ dnorm(0, 0.001)
  sigma_log_y ~ dunif(0,60)
  mu_Y <- exp(mu_log_y+sigma_log_y^2/2)
  sigma_Y <- sqrt(exp(2*mu_log_y+sigma_log_y^2)*(exp(sigma_log_y^2)-1))
  # Likelihood
  for(i in 1:n){
  y.SLA[i]~dlnorm(mu_log_y,1/sigma_log_y^2)}
  }
  ", fill=T)
  sink()
}

# Set intial conditions for the parameters within each chain 
inits = list(
list(
mu_log_y = 0,
sigma_log_y = .5
),
list(
mu_log_y = 1,
sigma_log_y = 1
)
)

# Write JAGS model and run

out.pooled.SLA = run.jags(model="Garbowski_Pooled_POPS_JAGS.R", method="rjparallel",
monitor=c( "mu_Y", "sigma_Y"),
data=data, n.chains=2, inits=inits,
adapt=250, burnin=250, sample=5000, thin=1)
output.pooled.SLA <- summary(out.pooled.SLA$mcmc)
# Trace plots chains for model parameters
plot(out.pooled.SLA$mcmc)
# Summary table 
as.table(output.pooled.SLA$statistics)
# Test for convergence using the Gelman diagnostic
gelman.diag(out.pooled.SLA$mcmc)

#
#
#
#
#
# Model with estimates of mu for each population - Measurements from ALL sample points are being lumped together 

y.POP.n = length(unique(HECO_data$POP_ID)) 

# Input Data 
data = list(
y.SLA = (HECO_data$SLA.TOT),
y.POP.group = HECO_data$POP_ID.index,
y.POP.n = length(unique(HECO_data$POP_ID)),
n = nrow(HECO_data))

# Initial conditions

inits = list(
list(
mu_pop = rep(5, y.POP.n),
mu_log = 0,
sigma_log = .5
),
list(
mu_pop = rep(10, y.POP.n),
mu_log = 1,
sigma_log = 1
)
)

{
  sink("Garbowski_Individual_POPS_JAGS.R")
  cat("
  model{
  # Priors
  

  mu_log ~ dnorm(0, 0.001)
  sigma_log ~ dunif(0,60)
  
  mu_y= exp(mu_log+sigma_log^2/2)
  sigma_y = sqrt(exp(2*mu_log+sigma_log^2)*(exp(sigma_log^2)-1))
  tau_y = 1/(sigma_y^2)
  
  # Likelihood
  for(i in 1:n){
  mu[i] = mu_pop[y.POP.group[i]]  
  y.SLA[i] ~ dnorm(mu[i], tau_y) }
  for (j in 1:y.POP.n){ #set loop to cycle through POPS
  mu_pop[j]~dlnorm(mu_y,tau_y)
  }}
  ", fill=T)
  sink()
}

# Write JAGS model and run
out.individual.SLA = run.jags(model="Garbowski_Individual_POPS_JAGS.R", method="rjparallel",
monitor=c("mu_pop"),
data=data, n.chains=2, inits=inits,
adapt=250, burnin=250, sample=5000, thin=1)
output.individual.SLA <- summary(out.individual.SLA$mcmc)
#  Trace plots of the chains for model parameters
plot(out.individual.SLA$mcmc)
# Summary table for the parameters
as.data.frame(output.individual.SLA$statistics)
# Test for convergence using the Gelman diagnostic
gelman.diag(out.individual.SLA$mcmc)

# Check against "regular" means
HECO_SLA<-HECO_data[, c("POP_ID","SLA.TOT" )]
aggregate(HECO_SLA[,2], list(HECO_SLA$POP_ID), mean)

#
#
#
#
# Model with estimates of mu for each population at 24 day timepoint. 
# Subset of data from the 24-day measurement being used.
# Estimates of variance by population included in this model.  
# It would be great to include data from all timepoints for a given species in the same model. 

HECO_24Day_data<-HECO_data [(HECO_data$Days%in%c("24")), ]
y.POP.n = length(unique(HECO_24Day_data$POP_ID)) 

# Input Data 
data = list(
  y.SLA = (HECO_24Day_data$SLA.TOT),
  y.POP.group = HECO_24Day_data$POP_ID.index,
  y.POP.n = length(unique(HECO_24Day_data$POP_ID)),
  n = nrow(HECO_24Day_data))

# Initial conditions

inits = list(
  list(
    mu_group = rep(5, y.POP.n),
    sigma_group = rep(10, y.POP.n),
    mu_log = 0,
    sigma_log = .5
  ),
  list(
    mu_group = rep(10, y.POP.n),
    sigma_group = rep(15, y.POP.n),
    mu_log = 1,
    sigma_log = 1
  )
)

{
  sink("Garbowski_Individual_POPS_JAGS.R")
  cat("
      model{
      # Priors
      
      mu_log ~ dnorm(0, 0.001)
      sigma_log ~ dunif(0,60)
      
      mu_y= exp(mu_log+sigma_log^2/2)
      sigma_y = sqrt(exp(2*mu_log+sigma_log^2)*(exp(sigma_log^2)-1))
      tau_y = 1/(sigma_y^2)
      
      # Likelihood
      for(i in 1:n){
      mu[i] = mu_group[y.POP.group[i]] 
      sigma[i] = sigma_group[y.POP.group[i]]
      tau [i] = 1/sigma [i]^2
      y.SLA[i] ~ dnorm(mu[i], tau [i]) }

      for (j in 1:y.POP.n){ # set loop to cycle through POPS
      mu_group[j] ~ dlnorm(mu_y, tau_y) # To make pops TOTALLY
      sigma_group[j] ~ dunif(0, 60) 
      }}
      ", fill=T)
  sink()
}

# Write JAGS model and run
out.individual.SLA = run.jags(model="Garbowski_Individual_POPS_JAGS.R", method="rjparallel",
                              monitor=c("mu_group", "sigma_group"),
                              data=data, n.chains=2, inits=inits,
                              adapt=250, burnin=250, sample=5000, thin=1)
output.individual.SLA <- summary(out.individual.SLA$mcmc)
#Produce trace plots of the chains for model parameters
plot(out.individual.SLA$mcmc)
#Produce a summary table for the parameters
as.data.frame(output.individual.SLA$statistics)
#Test for convergence using the Gelman diagnostic
gelman.diag(out.individual.SLA$mcmc)

#Check against "regular" means
HECO_SLA<-HECO_24Day_data[, c("POP_ID","SLA.TOT" )]
aggregate(HECO_SLA[,2], list(HECO_SLA$POP_ID), mean)