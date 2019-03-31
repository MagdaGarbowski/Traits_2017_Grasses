# Garbowski - 2017 CPNPP Traits
# Questions for Daniel Laughlin 
# March 30, 2019 

# Goals/What I'd appreciate help with (1 -3, maybe 4 down the road).#<br />
  #Technical Questions:#
  
# 1. How can I model variance (sigma or tau) from each population in JAGS? I have this figured out in STAN but can't seem to get it here. 
# 2. How can I include sampling date (harvest time, "H_num"" or "Days" in the data file) in my model
# I'd like to be able to compare estimates within/among population(s) but at different times? 
# Is the model below appropriate for all traits? Are the priors/distributions reasonable? Even for traits with values 0-1?
# Evenutally...How to compare among species?
  
# Examples of my process so far (translated to the best of my ability from STAN to JAGS).
# I am still "hacking" to try to get things to work 

# Load packages
library(ggplot2)
library(rjags)
library(parallel)
library(runjags)
library(reshape)
library(coda)

# Load data 

Trait2017_data<-read.csv("/Users/MagdaGarbowski/Traits_2017_Grasses/HECO_ELTR_PLPA_VUOC.csv")
# Subset Heterostipa comata (HECO) 
HECO_data<-Trait2017_data[Trait2017_data$SPECIES=="HECO",]
# Create integer index for POP_NUM
HECO_data$POP_ID.index<-as.integer(
factor(HECO_data[,c("POP_ID")],levels=c("HECO_AZNC_KV","HECO_UTC_KRR","HECO_UTC_TMR","HECO_UTEC_GR","HECO_UTSW_SP")))
HECO_data<-HECO_data[complete.cases(HECO_data), ]


# Function to link the sequential indices used in JAGS to the groups (POP_ID and H_num) in the data.

POP_ID_from_index = function(POP_ID, POP_ID.index, output ){
# POP_ID is a vector of populations (POP_ID)
# POP_ID.index is vector of sequential indices to POP_ID 

a = unique(as.vector(POP_ID)) 
b = unique(POP_ID.index)
POP_ID.key=as.data.frame(t(rbind(a,b))) # columns containing POP_ID indices paired with POP_ID number
names(POP_ID.key)= c(names(as.data.frame(POP_ID)), names(as.data.frame(POP_ID.index))) 
link.column.name = names(POP_ID.key)[2] # name of column for merging output data with POP_IDs
output2 = cbind(seq(1,nrow(output)),output) # give the output data sequential index the same as the POP_ID index
colnames(output2)[1]=link.column.name
POP_ID.data=as.data.frame(merge(POP_ID.key, output2, by = link.column.name )) # merge JAGS output with POP_ID names
return(POP_ID.data)
}


#  Pooled Model - Pooled data for all populations of HECO. Working with SLA column.

# Input Data
data = list(
y.SLA = (HECO_data$SLA.TOT),
n = length(HECO_data$SLA.TOT)) # number of observations (n) for the JAGS counter

# Model
{
  sink("Pooled_POPS_JAGS.R")
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

out.pooled.SLA = run.jags(model="Pooled_POPS_JAGS.R", method="rjparallel",
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


# Model with estimates of mu (intercepts) for each population. 
Measurements from ALL sample points are being lumped together 
```{r include=T, message=F, warning=F}
y.POP.n = length(unique(HECO_data$POP_ID)) # This is being difficult in the data call. Magda: fix this.
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
sigma = 10,
mu_log = 0,
sigma_log = .5
),
list(
mu_pop = rep(10, y.POP.n),
sigma = 15,
mu_log = 1,
sigma_log = 1
)
)

{
  sink("Individual_POPS_JAGS.R")
  cat("
  model{
  # Priors
  
  sigma ~ dunif (0,60)
  
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
out.individual.SLA = run.jags(model="Individual_POPS_JAGS.R", method="rjparallel",
monitor=c("mu_pop"),
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
HECO_SLA<-HECO_data[, c("POP_ID","SLA.TOT" )]
aggregate(HECO_SLA[,2], list(HECO_SLA$POP_ID), mean)

# Model with estimates of mu (intercepts) for each population at 24 day timepoint 
# I'm using a subset of data from the 24-day measurement. 
# It would be great to not do this and have all data from a species in the same model. 
# I assume I need to let the slopes vary by time (H_num). Add a time term?

HECO_24Day_data<-HECO_data [(HECO_data$Days%in%c("24")), ]

y.POP.n = length(unique(HECO_24Day_data$POP_ID)) # This is being difficult in the data = list. Magda: fix this.
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
    sigma = 10,
    mu_log = 0,
    sigma_log = .5
  ),
  list(
    mu_group = rep(10, y.POP.n),
    sigma = 15,
    mu_log = 1,
    sigma_log = 1
  )
)

{
  sink("Individual_POPS_JAGS.R")
  cat("
      model{
      # Priors
      
      sigma ~ dunif (0,60)
      
      mu_log ~ dnorm(0, 0.001)
      sigma_log ~ dunif(0,60)
      
      mu_y= exp(mu_log+sigma_log^2/2)
      sigma_y = sqrt(exp(2*mu_log+sigma_log^2)*(exp(sigma_log^2)-1))
      tau_y = 1/(sigma_y^2)
      
      # Likelihood
      for(i in 1:n){
      mu[i] = mu_group[y.POP.group[i]]  
      y.SLA[i] ~ dnorm(mu[i], tau_y) }
      for (j in 1:y.POP.n){ #set loop to cycle through POPS
      mu_group[j]~dlnorm(mu_y,tau_y)
      }}
      ", fill=T)
  sink()
}

# Write JAGS model and run
out.individual.SLA = run.jags(model="Individual_POPS_JAGS.R", method="rjparallel",
                              monitor=c("mu_group"),
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


# ATTEMPT at model with estimates of mu (intercepts) for each population varying with time

y.POP.n = length(unique(HECO_data$POP_ID)) # This is being difficult in the data call. Magda: fix this.
y.time.n = length(unique(HECO_data$H_num))
# Input Data 
data = list(
  y.SLA = (HECO_data$SLA.TOT),
  y.POP.group = HECO_data$POP_ID.index,
  y.POP.n = length(unique(HECO_data$POP_ID)),
  n = nrow(HECO_data))
  y.time = (HECO_data$H_num)
  y.time.n = length(unique(HECO_data$H_num))
  

# Initial conditions

inits = list(
  list(
    mu_pop = rep(5, y.POP.n),
    time = rep(.5,y.time.n),
    sigma = 10,
    mu_log = 0,
    sigma_log = .5
    eta = .2,
    kappa = .5
  ),
  list(
    mu_pop = rep(10, y.POP.n),
    time = rep(1,y.time.n),
    sigma = 15,
    mu_log = 1,
    sigma_log = 1
    eta = .2,
    kappa = 5
  )
)

{
  sink("Individual_Time_POPS_JAGS.R")
  cat("
      model{
      # Priors
      
      sigma ~ dunif (0,60)
      eta ~ dnorm (0, 1/1000^2) # change in mu? 
      kappa ~ dnorm (0, 1/1000^2) # when t = 0 what would mu be? Is this even relevant? 
      mu_time ~ dnorm(0, 1/1000^2)
      sigma_time ~ dunif (0, 200)
      
      mu_log ~ dnorm(0, 0.001)
      sigma_log ~ dunif(0,60)
      
      mu_y= exp(mu_log+sigma_log^2/2)
      sigma_y = sqrt(exp(2*mu_log+sigma_log^2)*(exp(sigma_log^2)-1))
      tau_y = 1/(sigma_y^2)
      tau_time = 1/sigma_time
      
      # Likelihood
      for(i in 1:n){
      mu[i] = mu_pop[y.POP.group[i]]  + time[y.time[i]]
      y.SLA[i] ~ dnorm(mu[i], tau_y) }

      # model for time effect (t's) (slope) and population (intercept) (j) on mu 
  for (t in 1:y.n.time){                     #set loop to calculate through time
      time [t] ~ dnorm (mu_time, tau_time)
}
  for (j in 1:y.POP.n){ #set loop to cycle through POPS
      mu_pop.t[j] = XXX    # average mu at time t for population j
      mu_group[j]~dlnorm(mu_y,tau_y)
      }}
      ", fill=T)
  sink()
}

# Write JAGS model and run
out.individual.SLA = run.jags(model="Individual_POPS_JAGS.R", method="rjparallel",
                              monitor=c("mu_pop"),
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
HECO_SLA<-HECO_data[, c("POP_ID","SLA.TOT" )]
aggregate(HECO_SLA[,2], list(HECO_SLA$POP_ID), mean)

