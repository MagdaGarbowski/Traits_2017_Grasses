### Traits analyses code
### 03-03-2019 
## Code below (with associated stan code) is for indidivudal species at individual time points. 
## Compares within a timepoint, among populations (2)

library(rstan)
library(plyr)
library(data.table)
library(ggplot2)
library(gridExtra)

SpeciesData<-read.csv("HECO_ELTR_PLPA_VUOC.csv")


#Total Weight = TotWD
#Total Root Length=SumOfLength.cm.
#SLA = SLA.Tot
#RMR = RMR
#SRL = S.RL
#LDMC = LDMC
#Height = HT

# Select out 10 days 
SpeciesData_10Day<-as.data.frame(SpeciesData [ which(SpeciesData$Days=='10'),  ])
SpeciesData_24Day<-as.data.frame(SpeciesData [ which(SpeciesData$Days=='24'),  ])
SpeciesData_42Day<-as.data.frame(SpeciesData [ which(SpeciesData$Days=='42'),  ])

####
####
####
####
####
####
### HECO

HECO.data<-as.data.frame(SpeciesData [ which(SpeciesData$SPECIES=='HECO'),  ])

### Duplicate POP_ID and change to PopNUM
HECO.data$PopNUM<-HECO.data$POP_ID
HECO.data$PopNUM<-as.integer(
  factor(HECO.data[,c("PopNUM")],levels=c("HECO_AZNC_KV","HECO_UTC_KRR","HECO_UTC_TMR","HECO_UTEC_GR","HECO_UTSW_SP")))

### 
HECO.RL = as.data.frame(HECO.data[,c(18,3,5)])
HECO.RL$PopNUM<-as.integer(HECO.RL$PopNUM)

### 10-Days 
HECO.RL_10<-HECO.RL [(HECO.RL$Days%in%c("10")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=4,
                 alpha_sigma=0,
                 beta_sigma=0.8,
                 N=nrow(HECO.RL_10),
                 n_grps = length(unique(HECO.RL_10$PopNUM)),
                 grp = HECO.RL_10$PopNUM,
                 y = HECO.RL_10$S.RL)
# Sample
HECO_fit_SRL_10day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SRL_10day, pars = c("mu", "sigma", "mu_regional"), digits = 2)
HECO_10day_SRL_mu_plot<-plot(HECO_fit_SRL_10day, pars = "mu",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_10day_SRL_sigma_plot<-plot(HECO_fit_SRL_10day, pars = "sigma",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_10day_SRL_sigma_plot<-plot(HECO_fit_SRL_10day, pars = "sigma",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))


# Check against empirical values
aggregate(HECO.RL_10[,3], by=list(PoppNUM=HECO.RL_10$PopNUM), FUN=mean)
aggregate(HECO.RL_10[,3], by=list(PoppNUM=HECO.RL_10$PopNUM), FUN=sd)

### 24-Days 
HECO.RL_24<-HECO.RL [(HECO.RL$Days%in%c("24")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.35,
                 N=nrow(HECO.RL_24),
                 n_grps = length(unique(HECO.RL_24$PopNUM)),
                 grp = HECO.RL_24$PopNUM,
                 y = HECO.RL_24$S.RL)
# Sample
HECO_fit_SRL_24day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SRL_24day, pars = c("mu", "sigma"))
HECO_24day_SRL_mu_plot<-plot(HECO_fit_SRL_24day, pars = "mu",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_24day_SRL_sigma_plot<-plot(HECO_fit_SRL_24day, pars = "sigma",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.RL_24[,3], by=list(PoppNUM=HECO.RL_24$PopNUM), FUN=mean)
aggregate(HECO.RL_24[,3], by=list(PoppNUM=HECO.RL_24$PopNUM), FUN=sd)

### 42-Days 
HECO.RL_42<-HECO.RL [(HECO.RL$Days%in%c("42")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.35,
                 N=nrow(HECO.RL_42),
                 n_grps = length(unique(HECO.RL_42$PopNUM)),
                 grp = HECO.RL_42$PopNUM,
                 y = HECO.RL_42$S.RL)
# Sample
HECO_fit_SRL_42day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SRL_42day, pars = c("mu", "sigma"))
HECO_42day_SRL_mu_plot<-plot(HECO_fit_SRL_42day, pars="mu",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_42day_SRL_sigma_plot<-plot(HECO_fit_SRL_42day, pars = "sigma",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.RL_42[,3], by=list(PoppNUM=HECO.RL_42$PopNUM), FUN=mean)
aggregate(HECO.RL_42[,3], by=list(PoppNUM=HECO.RL_42$PopNUM), FUN=sd)


### TOTAL ROOT LENGTH
HECO.SumOfLength.cm. = as.data.frame(HECO.data[,c(55,3,45)])
HECO.SumOfLength.cm.$PopNUM<-as.integer(HECO.SumOfLength.cm.$PopNUM)

### 10-Days 
HECO.SumOfLength.cm._10<-HECO.SumOfLength.cm. [(HECO.SumOfLength.cm.$Days%in%c("10")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,                  
                 alpha_sigma=0,                  
                 beta_sigma=.35,
                 N=nrow(HECO.SumOfLength.cm._10),
                 n_grps = length(unique(HECO.SumOfLength.cm._10$PopNUM)),
                 grp = HECO.SumOfLength.cm._10$PopNUM,
                 y = HECO.SumOfLength.cm._10$SumOfLength.cm.)
# Sample
HECO_fit_TRL_10day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_TRL_10day, pars = c("mu", "sigma"))
HECO_10day_TRL_mu_plot<-plot(HECO_fit_TRL_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_10day_TRL_sigma_plot<-plot(HECO_fit_TRL_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))


# Check against empirical values
aggregate(HECO.SumOfLength.cm._10[,2], by=list(PoppNUM=HECO.SumOfLength.cm._10$PopNUM), FUN=mean)
aggregate(HECO.SumOfLength.cm._10[,2], by=list(PoppNUM=HECO.SumOfLength.cm._10$PopNUM), FUN=sd)

### 24-Days 
HECO.SumOfLength.cm._24<-HECO.SumOfLength.cm. [(HECO.SumOfLength.cm.$Days%in%c("24")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=5,
                 alpha_sigma=0, 
                 beta_sigma=.2,
                 N=nrow(HECO.SumOfLength.cm._24),
                 n_grps = length(unique(HECO.SumOfLength.cm._24$PopNUM)),
                 grp = HECO.SumOfLength.cm._24$PopNUM,
                 y = HECO.SumOfLength.cm._24$SumOfLength.cm.)
# Sample
HECO_fit_TRL_24day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_TRL_24day, pars = c("mu", "sigma"))
HECO_24day_TRL_mu_plot<-plot(HECO_fit_TRL_24day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_24day_TRL_sigma_plot<-plot(HECO_fit_TRL_24day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.SumOfLength.cm._24[,2], by=list(PoppNUM=HECO.SumOfLength.cm._24$PopNUM), FUN=mean)
aggregate(HECO.SumOfLength.cm._24[,2], by=list(PoppNUM=HECO.SumOfLength.cm._24$PopNUM), FUN=sd)

### 42-Days 
HECO.SumOfLength.cm._42<-HECO.SumOfLength.cm. [(HECO.SumOfLength.cm.$Days%in%c("42")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0, 
                 beta_mu=20,                
                 alpha_sigma=0,               
                 beta_sigma=.25,
                 N=nrow(HECO.SumOfLength.cm._42),
                 n_grps = length(unique(HECO.SumOfLength.cm._42$PopNUM)),
                 grp = HECO.SumOfLength.cm._42$PopNUM,
                 y = HECO.SumOfLength.cm._42$SumOfLength.cm.)
# Sample
HECO_fit_TRL_42day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_TRL_42day, pars = c("mu", "sigma"))
HECO_42day_TRL_mu_plot<-plot(HECO_fit_TRL_42day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_42day_TRL_sigma_plot<-plot(HECO_fit_TRL_42day,outer_level=0.9, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.SumOfLength.cm._42[,2], by=list(PoppNUM=HECO.SumOfLength.cm._42$PopNUM), FUN=mean)
aggregate(HECO.SumOfLength.cm._42[,2], by=list(PoppNUM=HECO.SumOfLength.cm._42$PopNUM), FUN=sd)

### SPECIFIC LEAF AREA 
HECO.SLA. = as.data.frame(HECO.data[,c(55,48,45)])
HECO.SLA.$PopNUM<-as.integer(HECO.SLA.$PopNUM)

### 10-Days 
HECO.SLA._10<-HECO.SLA. [(HECO.SLA.$Days%in%c("10")), ]
HECO.SLA._10<-HECO.SLA._10[complete.cases(HECO.SLA._10), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=1.3,
                 alpha_sigma=0,
                 beta_sigma=.3,
                 N=nrow(HECO.SLA._10),
                 n_grps = length(unique(HECO.SLA._10$PopNUM)),
                 grp = HECO.SLA._10$PopNUM,
                 y = HECO.SLA._10$SLA.TOT)

# Sample
HECO_fit_SLA_10day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SLA_10day, pars = c("mu", "sigma"))
HECO_10day_SLA_mu_plot<-plot(HECO_fit_SLA_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_10day_SLA_sigma_plot<-plot(HECO_fit_SLA_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.SLA._10[,2], by=list(PoppNUM=HECO.SLA._10$PopNUM), FUN=mean)
aggregate(HECO.SLA._10[,2], by=list(PoppNUM=HECO.SLA._10$PopNUM), FUN=sd)


### 24-Days 
HECO.SLA._24<-HECO.SLA. [(HECO.SLA.$Days%in%c("24")), ]
HECO.SLA._24<-HECO.SLA._24[complete.cases(HECO.SLA._24), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=1.8,     
                 alpha_sigma=0,     
                 beta_sigma=.2,
                 N=nrow(HECO.SLA._24),
                 n_grps = length(unique(HECO.SLA._24$PopNUM)),
                 grp = HECO.SLA._24$PopNUM,
                 y = HECO.SLA._24$SLA.TOT)
# Sample
HECO_fit_SLA_24day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SLA_24day, pars = c("mu", "sigma"))
HECO_24day_SLA_mu_plot<-plot(HECO_fit_SLA_24day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_24day_SLA_sigma_plot<-plot(HECO_fit_SLA_24day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.SLA._24[,2], by=list(PoppNUM=HECO.SLA._24$PopNUM), FUN=mean)
aggregate(HECO.SLA._24[,2], by=list(PoppNUM=HECO.SLA._24$PopNUM), FUN=sd)

### 42-Days 
HECO.SLA._42<-HECO.SLA. [(HECO.SLA.$Days%in%c("42")), ]
HECO.SLA._42<-HECO.SLA._42[complete.cases(HECO.SLA._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.22,
                 N=nrow(HECO.SLA._42),
                 n_grps = length(unique(HECO.SLA._42$PopNUM)),
                 grp = HECO.SLA._42$PopNUM,
                 y = HECO.SLA._42$SLA.TOT)
# Sample
HECO_fit_SLA_42day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SLA_42day, pars = c("mu", "sigma"))
HECO_42day_SLA_mu_plot<-plot(HECO_fit_SLA_42day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_42day_SLA_sigma_plot<-plot(HECO_fit_SLA_42day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.SLA._42[,2], by=list(PoppNUM=HECO.SLA._42$PopNUM), FUN=mean)
aggregate(HECO.SLA._42[,2], by=list(PoppNUM=HECO.SLA._42$PopNUM), FUN=sd)

### HEIGHT
HECO.HT. = as.data.frame(HECO.data[,c(55,19,45)])
HECO.HT.$PopNUM<-as.integer(HECO.HT.$PopNUM)

### 10-Days 
HECO.HT._10<-HECO.HT. [(HECO.HT.$Days%in%c("10")), ]
HECO.HT._10<-HECO.HT._10[complete.cases(HECO.HT._10), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.12,
                 N=nrow(HECO.HT._10),
                 n_grps = length(unique(HECO.HT._10$PopNUM)),
                 grp = HECO.HT._10$PopNUM,
                 y = HECO.HT._10$HT)

# Sample
HECO_fit_HT_10day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_HT_10day, pars = c("mu", "sigma"))
HECO_10day_HT_mu_plot<-plot(HECO_fit_HT_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_10day_HT_sigma_plot<-plot(HECO_fit_HT_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.HT._10[,2], by=list(PoppNUM=HECO.HT._10$PopNUM), FUN=mean)
aggregate(HECO.HT._10[,2], by=list(PoppNUM=HECO.HT._10$PopNUM), FUN=sd)

### 24-Days 
HECO.HT._24<-HECO.HT. [(HECO.HT.$Days%in%c("24")), ]
HECO.HT._24<-HECO.HT._24[complete.cases(HECO.HT._24), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.5,
                 alpha_sigma=0,
                 beta_sigma=.2,
                 N=nrow(HECO.HT._24),
                 n_grps = length(unique(HECO.HT._24$PopNUM)),
                 grp = HECO.HT._24$PopNUM,
                 y = HECO.HT._24$HT)
# Sample
HECO_fit_HT_24day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_HT_24day, pars = c("mu", "sigma"))
HECO_24day_HT_mu_plot<-plot(HECO_fit_HT_24day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_24day_HT_sigma_plot<-plot(HECO_fit_HT_24day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.HT._24[,2], by=list(PoppNUM=HECO.HT._24$PopNUM), FUN=mean)
aggregate(HECO.HT._24[,2], by=list(PoppNUM=HECO.HT._24$PopNUM), FUN=sd)

### 42-Days 
HECO.HT._42<-HECO.HT. [(HECO.HT.$Days%in%c("42")), ]
HECO.HT._42<-HECO.HT._42[complete.cases(HECO.HT._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=4,  
                 alpha_sigma=0,
                 beta_sigma=.2,
                 N=nrow(HECO.HT._42),
                 n_grps = length(unique(HECO.HT._42$PopNUM)),
                 grp = HECO.HT._42$PopNUM,
                 y = HECO.HT._42$HT)
# Sample
HECO_fit_HT_42day = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_HT_42day, pars = c("mu", "sigma"))
HECO_42day_HT_mu_plot<-plot(HECO_fit_HT_42day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))
HECO_42day_HT_sigma_plot<-plot(HECO_fit_HT_42day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "forestgreen", "deepskyblue3"))

# Check against empirical values
aggregate(HECO.HT._42[,2], by=list(PoppNUM=HECO.HT._42$PopNUM), FUN=mean)
aggregate(HECO.HT._42[,2], by=list(PoppNUM=HECO.HT._42$PopNUM), FUN=sd)

```

##### H. comata - Total root length (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(HECO_10day_TRL_mu_plot+xlim(10,200), HECO_10day_TRL_sigma_plot+xlim(0,200),
             HECO_24day_TRL_mu_plot+xlim(10,500), HECO_24day_TRL_sigma_plot+xlim(0,500),
             HECO_42day_TRL_mu_plot+xlim(10,500), HECO_42day_TRL_sigma_plot+xlim(0,500),
             nrow=3 )

```

##### H. comata - Specific root length (cm/g)
```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(HECO_10day_SRL_mu_plot+xlim(10,380), HECO_10day_SRL_sigma_plot+xlim(0,300),
             HECO_24day_SRL_mu_plot+xlim(10,380), HECO_24day_SRL_sigma_plot+xlim(0,300),
             HECO_42day_SRL_mu_plot+xlim(10,380), HECO_42day_SRL_sigma_plot+xlim(0,300),
             nrow=3 )

```

##### H. comata - Specific leaf area (cm^2/g)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(HECO_10day_SLA_mu_plot+xlim(0,20), HECO_10day_SLA_sigma_plot+xlim(0,10),
             HECO_24day_SLA_mu_plot+xlim(0,20), HECO_24day_SLA_sigma_plot+xlim(0,10),
             HECO_42day_SLA_mu_plot+xlim(0,50), HECO_42day_SLA_sigma_plot+xlim(0,30),
             nrow=3 )

```

##### H. comata - Height (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(HECO_10day_HT_mu_plot+xlim(0,35), HECO_10day_HT_sigma_plot+xlim(0,20),
             HECO_24day_HT_mu_plot+xlim(0,35), HECO_24day_HT_sigma_plot+xlim(0,20),
             HECO_42day_HT_mu_plot+xlim(0,35), HECO_42day_HT_sigma_plot+xlim(0,20),
             nrow=3 )

```

##### "Regional" Hesperostipa comata SLA at 24 days and 42 days 
```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, results='hide'}
### SPECIFIC LEAF AREA 
HECO.SLA. = as.data.frame(HECO.data[,c(55,48,45)])
HECO.SLA.$PopNUM<-as.integer(HECO.SLA.$PopNUM)

### 24-Days 
HECO.SLA._24<-HECO.SLA. [(HECO.SLA.$Days%in%c("24")), ]
HECO.SLA._24<-HECO.SLA._24[complete.cases(HECO.SLA._24), ]
### Add a PopALL Column = 1 to change groups to PopALL
HECO.SLA._24$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(HECO.SLA._24),
                 n_grps = length(unique(HECO.SLA._24$PopALL)),
                 grp = HECO.SLA._24$PopALL,
                 y = HECO.SLA._24$SLA.TOT)

# Sample
HECO_fit_SLA_24day_PopALL = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SLA_24day_PopALL, pars = c("mu", "sigma"))
HECO_24day_SLA_mu_plot_PopALL<-plot(HECO_fit_SLA_24day_PopALL, pars = "mu", fill_color="black")
HECO_24day_SLA_sigma_plot_PopALL<-plot(HECO_fit_SLA_24day_PopALL, pars = "sigma", fill_color=c("black"))
HECO_24day_SLA_mu_PopALL<-plot(HECO_fit_SLA_24day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(HECO.SLA._24[,2], by=list(PoppNUM=HECO.SLA._24$PopALL), FUN=mean)
aggregate(HECO.SLA._24[,2], by=list(PoppNUM=HECO.SLA._24$PopALL), FUN=sd)




### 42-Days 
HECO.SLA._42<-HECO.SLA. [(HECO.SLA.$Days%in%c("42")), ]
HECO.SLA._42<-HECO.SLA._42[complete.cases(HECO.SLA._42), ]
### Add a PopALL Column = 1 to change groups to PopALL
HECO.SLA._42$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(HECO.SLA._42),
                 n_grps = length(unique(HECO.SLA._42$PopALL)),
                 grp = HECO.SLA._42$PopALL,
                 y = HECO.SLA._42$SLA.TOT)

# Sample
HECO_fit_SLA_42day_PopALL = sampling(mod, data = data_list)

# Look at output
print(HECO_fit_SLA_42day_PopALL, pars = c("mu", "sigma"))
HECO_42day_SLA_mu_plot_PopALL<-plot(HECO_fit_SLA_42day_PopALL, pars = "mu", fill_color="black")
HECO_42day_SLA_sigma_plot_PopALL<-plot(HECO_fit_SLA_42day_PopALL, pars = "sigma", fill_color=c("black"))
HECO_42day_SLA_mu_PopALL<-plot(HECO_fit_SLA_42day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(HECO.SLA._42[,2], by=list(PoppNUM=HECO.SLA._42$PopALL), FUN=mean)
aggregate(HECO.SLA._42[,2], by=list(PoppNUM=HECO.SLA._42$PopALL), FUN=sd)

```
```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE}

HECO_24day_SLA_mu_PopALL+geom_vline(xintercept = 9.72, linetype="solid", size=.7, color="black")+
  xlim(6,25)+
  geom_vline(xintercept = 15.34, linetype="dashed", size=.9, color="deepskyblue2")+
  geom_vline(xintercept = 9.83, linetype="dashed", size=.7, color="olivedrab3")+
  geom_vline(xintercept = 7.76, linetype="dashed", size=.7, color="olivedrab4")+
  geom_vline(xintercept = 9.08, linetype="dashed", size=.7, color="forestgreen")+
  geom_vline(xintercept = 9.75, linetype="dashed", size=.7, color="deepskyblue3")+
  geom_vline(xintercept = 23.75, linetype="dotted", size=.7, color="grey40")

HECO_42day_SLA_mu_PopALL+geom_vline(xintercept = 16.91, linetype="solid", size=.7, color="black")+
  xlim(6,25)+
  geom_vline(xintercept = 10.67, linetype="dashed", size=.9, color="deepskyblue2")+
  geom_vline(xintercept = 17.74, linetype="dashed", size=.7, color="olivedrab3")+
  geom_vline(xintercept = 17.35, linetype="dashed", size=.7, color="olivedrab4")+
  geom_vline(xintercept = 24.6, linetype="dashed", size=.7, color="forestgreen")+
  geom_vline(xintercept = 19.75, linetype="dashed", size=.7, color="deepskyblue3")+
  geom_vline(xintercept = 23.75, linetype="dotted", size=.7, color="grey40")


```


```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, results='hide'}
#
#
#
#
#
#
#
#
#
##
#
####
####
####
####
####
####
### PLPA

HEVI.data<-as.data.frame(SpeciesData [ which(SpeciesData$SPECIES=='HEVI'),  ])
### Duplicate POP_ID and change to PopNUM
HEVI.data$PopNUM<-HEVI.data$POP_ID

HEVI.data[,c("PopNUM")][HEVI.data[,c("PopNUM")]=="HEVI_AZNC_KV"]<-as.integer(1)
HEVI.data[,c("PopNUM")][HEVI.data[,c("PopNUM")]=="HEVI_COSW_DE"]<-as.integer(2)
HEVI.data[,c("PopNUM")][HEVI.data[,c("PopNUM")]=="HEVI_NMNW_NC"]<-as.integer(3)
HEVI.data[,c("PopNUM")][HEVI.data[,c("PopNUM")]=="HEVI_NMNW_ND"]<-as.integer(4)


HEVI.RL = as.data.frame(HEVI.data[,c(55,45,49)])
HEVI.RL$PopNUM<-as.integer(HEVI.RL$PopNUM)

### 10-Days 
#HEVI.RL_10<-HEVI.RL [(HEVI.RL$Days%in%c("10")), ]
#HEVI.RL_10<-HEVI.RL_10[complete.cases(HEVI.RL_10), ]
# Compile model
#mod = stan_model("SRL_10Day.stan")
#data_list = list(alpha_mu=0,
#beta_mu=2,                  
#alpha_sigma=0,                 
#beta_sigma=.35,
#N=nrow(HEVI.RL_10),
#  n_grps = length(unique(HEVI.RL_10$PopNUM)),
#   grp = HEVI.RL_10$PopNUM,
#  y = HEVI.RL_10$S.RL)
# Sample
#HEVI_fit_SRL_10day = sampling(mod, data = data_list)

# Look at output
#print(HEVI_fit_SRL_10day, pars = c("mu", "sigma"))
#plot(HEVI_fit_SRL_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
#plot(HEVI_fit_SRL_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
#aggregate(HEVI.RL_10[,3], by=list(PoppNUM=HEVI.RL_10$PopNUM), FUN=mean)
#aggregate(HEVI.RL_10[,3], by=list(PoppNUM=HEVI.RL_10$PopNUM), FUN=sd)

### 24-Days 
HEVI.RL_24<-HEVI.RL [(HEVI.RL$Days%in%c("24")), ]
HEVI.RL_24<-HEVI.RL_24[complete.cases(HEVI.RL_24), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.35,
                 N=nrow(HEVI.RL_24),
                 n_grps = length(unique(HEVI.RL_24$PopNUM)),
                 grp = HEVI.RL_24$PopNUM,
                 y = HEVI.RL_24$S.RL)
# Sample
HEVI_fit_SRL_24day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_SRL_24day, pars = c("mu", "sigma"))
HEVI_24day_SRL_mu_plot<-plot(HEVI_fit_SRL_24day, pars = "mu",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_24day_SRL_sigma_plot<-plot(HEVI_fit_SRL_24day, pars = "sigma",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.RL_24[,3], by=list(PoppNUM=HEVI.RL_24$PopNUM), FUN=mean)
aggregate(HEVI.RL_24[,3], by=list(PoppNUM=HEVI.RL_24$PopNUM), FUN=sd)

### 42-Days 
HEVI.RL_42<-HEVI.RL [(HEVI.RL$Days%in%c("42")), ]
HEVI.RL_42<-HEVI.RL_42[complete.cases(HEVI.RL_42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=1.8,
                 alpha_sigma=0,
                 beta_sigma=.35,
                 N=nrow(HEVI.RL_42),
                 n_grps = length(unique(HEVI.RL_42$PopNUM)),
                 grp = HEVI.RL_42$PopNUM,
                 y = HEVI.RL_42$S.RL)
# Sample
HEVI_fit_SRL_42day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_SRL_42day, pars = c("mu", "sigma"))
HEVI_42day_SRL_mu_plot<-plot(HEVI_fit_SRL_42day, pars="mu",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_42day_SRL_sigma_plot<-plot(HEVI_fit_SRL_42day, pars = "sigma",fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.RL_42[,3], by=list(PoppNUM=HEVI.RL_42$PopNUM), FUN=mean)
aggregate(HEVI.RL_42[,3], by=list(PoppNUM=HEVI.RL_42$PopNUM), FUN=sd)


### TOTAL ROOT LENGTH
HEVI.SumOfLength.cm. = as.data.frame(HEVI.data[,c(55,3,45)])
HEVI.SumOfLength.cm.$PopNUM<-as.integer(HEVI.SumOfLength.cm.$PopNUM)

### 10-Days 
HEVI.SumOfLength.cm._10<-HEVI.SumOfLength.cm. [(HEVI.SumOfLength.cm.$Days%in%c("10")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=3,                  
                 alpha_sigma=0,                  
                 beta_sigma=.25,
                 N=nrow(HEVI.SumOfLength.cm._10),
                 n_grps = length(unique(HEVI.SumOfLength.cm._10$PopNUM)),
                 grp = HEVI.SumOfLength.cm._10$PopNUM,
                 y = HEVI.SumOfLength.cm._10$SumOfLength.cm.)
# Sample
HEVI_fit_TRL_10day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_TRL_10day, pars = c("mu", "sigma"))
HEVI_10day_TRL_mu_plot<-plot(HEVI_fit_TRL_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_10day_TRL_sigma_plot<-plot(HEVI_fit_TRL_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))


# Check against empirical values
aggregate(HEVI.SumOfLength.cm._10[,2], by=list(PoppNUM=HEVI.SumOfLength.cm._10$PopNUM), FUN=mean)
aggregate(HEVI.SumOfLength.cm._10[,2], by=list(PoppNUM=HEVI.SumOfLength.cm._10$PopNUM), FUN=sd)

### 24-Days 
HEVI.SumOfLength.cm._24<-HEVI.SumOfLength.cm. [(HEVI.SumOfLength.cm.$Days%in%c("24")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0, 
                 beta_sigma=.3,
                 N=nrow(HEVI.SumOfLength.cm._24),
                 n_grps = length(unique(HEVI.SumOfLength.cm._24$PopNUM)),
                 grp = HEVI.SumOfLength.cm._24$PopNUM,
                 y = HEVI.SumOfLength.cm._24$SumOfLength.cm.)
# Sample
HEVI_fit_TRL_24day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_TRL_24day, pars = c("mu", "sigma"))
HEVI_24day_TRL_mu_plot<-plot(HEVI_fit_TRL_24day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_24day_TRL_sigma_plot<-plot(HEVI_fit_TRL_24day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.SumOfLength.cm._24[,2], by=list(PoppNUM=HEVI.SumOfLength.cm._24$PopNUM), FUN=mean)
aggregate(HEVI.SumOfLength.cm._24[,2], by=list(PoppNUM=HEVI.SumOfLength.cm._24$PopNUM), FUN=sd)

### 42-Days 
HEVI.SumOfLength.cm._42<-HEVI.SumOfLength.cm. [(HEVI.SumOfLength.cm.$Days%in%c("42")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0, 
                 beta_mu=4.2,                
                 alpha_sigma=0,               
                 beta_sigma=.3,
                 N=nrow(HEVI.SumOfLength.cm._42),
                 n_grps = length(unique(HEVI.SumOfLength.cm._42$PopNUM)),
                 grp = HEVI.SumOfLength.cm._42$PopNUM,
                 y = HEVI.SumOfLength.cm._42$SumOfLength.cm.)
# Sample
HEVI_fit_TRL_42day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_TRL_42day, pars = c("mu", "sigma"))
HEVI_42day_TRL_mu_plot<-plot(HEVI_fit_TRL_42day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_42day_TRL_sigma_plot<-plot(HEVI_fit_TRL_42day,outer_level=0.9, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.SumOfLength.cm._42[,2], by=list(PoppNUM=HEVI.SumOfLength.cm._42$PopNUM), FUN=mean)
aggregate(HEVI.SumOfLength.cm._42[,2], by=list(PoppNUM=HEVI.SumOfLength.cm._42$PopNUM), FUN=sd)

### SPECIFIC LEAF AREA 
HEVI.SLA. = as.data.frame(HEVI.data[,c(55,48,45)])
HEVI.SLA.$PopNUM<-as.integer(HEVI.SLA.$PopNUM)

### 10-Days 
#HEVI.SLA._10<-HEVI.SLA. [(HEVI.SLA.$Days%in%c("10")), ]
#HEVI.SLA._10<-HEVI.SLA._10[complete.cases(HEVI.SLA._10), ]

# Compile model
#mod = stan_model("Traits_1.stan")
#data_list = list(alpha_mu=0,
#  beta_mu=1.3,
#  alpha_sigma=0,
#  beta_sigma=.3,
#  N=nrow(HEVI.SLA._10),
#  n_grps = length(unique(HEVI.SLA._10$PopNUM)),
#  grp = HEVI.SLA._10$PopNUM,
#  y = HEVI.SLA._10$SLA.TOT)

# Sample
#HEVI_fit_SLA_10day = sampling(mod, data = data_list)

# Look at output
#print(HEVI_fit_SLA_10day, pars = c("mu", "sigma"))
#HEVI_10day_SLA_mu_plot<-plot(HEVI_fit_SLA_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
#HEVI_10day_SLA_sigma_plot<-plot(HEVI_fit_SLA_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
#aggregate(HEVI.SLA._10[,2], by=list(PoppNUM=HEVI.SLA._10$PopNUM), FUN=mean)
#aggregate(HEVI.SLA._10[,2], by=list(PoppNUM=HEVI.SLA._10$PopNUM), FUN=sd)


### 24-Days 
HEVI.SLA._24<-HEVI.SLA. [(HEVI.SLA.$Days%in%c("24")), ]
HEVI.SLA._24<-HEVI.SLA._24[complete.cases(HEVI.SLA._24), ]
### Remove popNUM = 2 because it only have one value 
#HEVI.SLA._24<-HEVI.SLA._24[!(HEVI.SLA._24$PopNUM%in%c("2")),]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=6,     
                 alpha_sigma=0,     
                 beta_sigma=.3,
                 N=nrow(HEVI.SLA._24),
                 n_grps = length(unique(HEVI.SLA._24$PopNUM)),
                 grp = HEVI.SLA._24$PopNUM,
                 y = HEVI.SLA._24$SLA.TOT)
# Sample
HEVI_fit_SLA_24day = sampling(mod, data = data_list)


# Look at output
print(HEVI_fit_SLA_24day, pars = c("mu", "sigma"))
HEVI_24day_SLA_mu_plot<-plot(HEVI_fit_SLA_24day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_24day_SLA_sigma_plot<-plot(HEVI_fit_SLA_24day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.SLA._24[,2], by=list(PoppNUM=HEVI.SLA._24$PopNUM), FUN=mean)
aggregate(HEVI.SLA._24[,2], by=list(PoppNUM=HEVI.SLA._24$PopNUM), FUN=sd)

### 42-Days 
HEVI.SLA._42<-HEVI.SLA. [(HEVI.SLA.$Days%in%c("42")), ]
HEVI.SLA._42<-HEVI.SLA._42[complete.cases(HEVI.SLA._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.25,
                 N=nrow(HEVI.SLA._42),
                 n_grps = length(unique(HEVI.SLA._42$PopNUM)),
                 grp = HEVI.SLA._42$PopNUM,
                 y = HEVI.SLA._42$SLA.TOT)
# Sample
HEVI_fit_SLA_42day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_SLA_42day, pars = c("mu", "sigma"))
HEVI_42day_SLA_mu_plot<-plot(HEVI_fit_SLA_42day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_42day_SLA_sigma_plot<-plot(HEVI_fit_SLA_42day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.SLA._42[,2], by=list(PoppNUM=HEVI.SLA._42$PopNUM), FUN=mean)
aggregate(HEVI.SLA._42[,2], by=list(PoppNUM=HEVI.SLA._42$PopNUM), FUN=sd)

### HEIGHT
HEVI.HT. = as.data.frame(HEVI.data[,c(55,19,45)])
HEVI.HT.$PopNUM<-as.integer(HEVI.HT.$PopNUM)

### 10-Days 
HEVI.HT._10<-HEVI.HT. [(HEVI.HT.$Days%in%c("10")), ]
HEVI.HT._10<-HEVI.HT._10[complete.cases(HEVI.HT._10), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.2,
                 N=nrow(HEVI.HT._10),
                 n_grps = length(unique(HEVI.HT._10$PopNUM)),
                 grp = HEVI.HT._10$PopNUM,
                 y = HEVI.HT._10$HT)

# Sample
HEVI_fit_HT_10day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_HT_10day, pars = c("mu", "sigma"))
HEVI_10day_HT_mu_plot<-plot(HEVI_fit_HT_10day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_10day_HT_sigma_plot<-plot(HEVI_fit_HT_10day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.HT._10[,2], by=list(PoppNUM=HEVI.HT._10$PopNUM), FUN=mean)
aggregate(HEVI.HT._10[,2], by=list(PoppNUM=HEVI.HT._10$PopNUM), FUN=sd)

### 24-Days 
HEVI.HT._24<-HEVI.HT. [(HEVI.HT.$Days%in%c("24")), ]
HEVI.HT._24<-HEVI.HT._24[complete.cases(HEVI.HT._24), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.5,
                 alpha_sigma=0,
                 beta_sigma=.15,
                 N=nrow(HEVI.HT._24),
                 n_grps = length(unique(HEVI.HT._24$PopNUM)),
                 grp = HEVI.HT._24$PopNUM,
                 y = HEVI.HT._24$HT)
# Sample
HEVI_fit_HT_24day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_HT_24day, pars = c("mu", "sigma"))
HEVI_24day_HT_mu_plot<-plot(HEVI_fit_HT_24day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_24day_HT_sigma_plot<-plot(HEVI_fit_HT_24day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.HT._24[,2], by=list(PoppNUM=HEVI.HT._24$PopNUM), FUN=mean)
aggregate(HEVI.HT._24[,2], by=list(PoppNUM=HEVI.HT._24$PopNUM), FUN=sd)

### 42-Days 
HEVI.HT._42<-HEVI.HT. [(HEVI.HT.$Days%in%c("42")), ]
HEVI.HT._42<-HEVI.HT._42[complete.cases(HEVI.HT._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=3,  
                 alpha_sigma=0,
                 beta_sigma=.15,
                 N=nrow(HEVI.HT._42),
                 n_grps = length(unique(HEVI.HT._42$PopNUM)),
                 grp = HEVI.HT._42$PopNUM,
                 y = HEVI.HT._42$HT)
# Sample
HEVI_fit_HT_42day = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_HT_42day, pars = c("mu", "sigma"))
HEVI_42day_HT_mu_plot<-plot(HEVI_fit_HT_42day, pars = "mu", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))
HEVI_42day_HT_sigma_plot<-plot(HEVI_fit_HT_42day, pars = "sigma", fill_color=c("deepskyblue2","olivedrab3", "olivedrab4", "olivedrab2"))

# Check against empirical values
aggregate(HEVI.HT._42[,2], by=list(PoppNUM=HEVI.HT._42$PopNUM), FUN=mean)
aggregate(HEVI.HT._42[,2], by=list(PoppNUM=HEVI.HT._42$PopNUM), FUN=sd)




```

##### H. villosa - Total root length (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(HEVI_10day_TRL_mu_plot+xlim(0,15), HEVI_10day_TRL_sigma_plot+xlim(0,15),
             HEVI_24day_TRL_mu_plot+xlim(0,60), HEVI_24day_TRL_sigma_plot+xlim(0,60),
             HEVI_42day_TRL_mu_plot+xlim(0,120), HEVI_42day_TRL_sigma_plot+xlim(0,120),
             nrow=3 )

```

##### H. villosa - Specific root length (cm/g)
```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
blankplot<-ggplot()
grid.arrange(blankplot,blankplot,HEVI_24day_SRL_mu_plot+xlim(10,380), HEVI_24day_SRL_sigma_plot+xlim(0,300),
             HEVI_42day_SRL_mu_plot+xlim(10,380), HEVI_42day_SRL_sigma_plot+xlim(0,300),
             nrow=3 )

```

##### H. villosa - Specific leaf area (cm^2/g)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
blankplot<-ggplot()
grid.arrange(blankplot,blankplot,
             HEVI_24day_SLA_mu_plot+xlim(0,60), HEVI_24day_SLA_sigma_plot+xlim(0,60),
             HEVI_42day_SLA_mu_plot+xlim(0,60), HEVI_42day_SLA_sigma_plot+xlim(0,60),
             nrow=3 )

```

##### H. villosa - Height (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(HEVI_10day_HT_mu_plot+xlim(0,4), HEVI_10day_HT_sigma_plot+xlim(0,4),
             HEVI_24day_HT_mu_plot+xlim(0,4), HEVI_24day_HT_sigma_plot+xlim(0,4),
             HEVI_42day_HT_mu_plot+xlim(0,4), HEVI_42day_HT_sigma_plot+xlim(0,4),
             nrow=3 )

```

##### H. villosa "regional" SLA at 24 days and 42 days
```{r, include=FALSE, message=FALSE, warning=FALSE, error=FALSE, echo=FALSE ,results='hide'}
### SPECIFIC LEAF AREA 
HEVI.SLA. = as.data.frame(HEVI.data[,c(55,48,45)])
HEVI.SLA.$PopNUM<-as.integer(HEVI.SLA.$PopNUM)

### 24-Days 
HEVI.SLA._24<-HEVI.SLA. [(HEVI.SLA.$Days%in%c("24")), ]
HEVI.SLA._24<-HEVI.SLA._24[complete.cases(HEVI.SLA._24), ]
### Add a PopALL Column = 1 to change groups to PopALL
HEVI.SLA._24$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(HEVI.SLA._24),
                 n_grps = length(unique(HEVI.SLA._24$PopALL)),
                 grp = HEVI.SLA._24$PopALL,
                 y = HEVI.SLA._24$SLA.TOT)

# Sample
HEVI_fit_SLA_24day_PopALL = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_SLA_24day_PopALL, pars = c("mu", "sigma"))
HEVI_24day_SLA_mu_plot_PopALL<-plot(HEVI_fit_SLA_24day_PopALL, pars = "mu", fill_color="black")
HEVI_24day_SLA_sigma_plot_PopALL<-plot(HEVI_fit_SLA_24day_PopALL, pars = "sigma", fill_color=c("black"))
HEVI_24day_SLA_mu_PopALL<-plot(HEVI_fit_SLA_24day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(HEVI.SLA._24[,2], by=list(PopNUM=HEVI.SLA._24$PopALL), FUN=mean)
aggregate(HEVI.SLA._24[,2], by=list(PopNUM=HEVI.SLA._24$PopALL), FUN=sd)



### 42-Days 
HEVI.SLA._42<-HEVI.SLA. [(HEVI.SLA.$Days%in%c("42")), ]

HEVI.SLA._42<-HEVI.SLA._42[complete.cases(HEVI.SLA._42), ]
### Add a PopALL Column = 1 to change groups to PopALL
HEVI.SLA._42$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(HEVI.SLA._42),
                 n_grps = length(unique(HEVI.SLA._42$PopALL)),
                 grp = HEVI.SLA._42$PopALL,
                 y = HEVI.SLA._42$SLA.TOT)

# Sample
HEVI_fit_SLA_42day_PopALL = sampling(mod, data = data_list)

# Look at output
print(HEVI_fit_SLA_42day_PopALL, pars = c("mu", "sigma"))
HEVI_42day_SLA_mu_plot_PopALL<-plot(HEVI_fit_SLA_42day_PopALL, pars = "mu", fill_color="black")
HEVI_42day_SLA_sigma_plot_PopALL<-plot(HEVI_fit_SLA_42day_PopALL, pars = "sigma", fill_color=c("black"))
HEVI_42day_SLA_mu_PopALL<-plot(HEVI_fit_SLA_42day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(HEVI.SLA._42[,2], by=list(PoppNUM=HEVI.SLA._42$PopALL), FUN=mean)
aggregate(HEVI.SLA._42[,2], by=list(PoppNUM=HEVI.SLA._42$PopALL), FUN=sd)
```

```{r, echo=FALSE, fig.width=6.6, fig.height=3.8, warning=FALSE, error=FALSE, message=FALSE}

HEVI_24day_SLA_mu_PopALL+geom_vline(xintercept = 10.5, linetype="solid", size=.7, color="black")+
  xlim(0,30)+
  geom_vline(xintercept = 11.34, linetype="dashed", size=.9, color="deepskyblue2")+
  geom_vline(xintercept = 4.4, linetype="dashed", size=.7, color="olivedrab3")+
  geom_vline(xintercept = 14.49, linetype="dashed", size=.7, color="olivedrab4")+
  geom_vline(xintercept = 19.64, linetype="dashed", size=.7, color="olivedrab2")


HEVI_42day_SLA_mu_PopALL+geom_vline(xintercept = 28.5, linetype="solid", size=.7, color="black")+
  xlim(0,45)+
  geom_vline(xintercept = 30.77, linetype="dashed", size=.9, color="deepskyblue2")+
  geom_vline(xintercept = 25.98, linetype="dashed", size=.7, color="olivedrab3")+
  geom_vline(xintercept = 39.91, linetype="dashed", size=.7, color="olivedrab4")+
  geom_vline(xintercept = 26.75, linetype="dashed", size=.7, color="olivedrab2")

```

```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, results='hide'}
#
#
#
#
#
#
#
#
#
##
#
####
####
####
####
####
####
### PLPA

PLPA.data<-as.data.frame(SpeciesData [ which(SpeciesData$SPECIES=='PLPA'),  ])
### Duplicate POP_ID and change to PopNUM
PLPA.data$PopNUM<-PLPA.data$POP_ID

PLPA.data[,c("PopNUM")][PLPA.data[,c("PopNUM")]=="PLPA_AZC_TC_"]<-as.integer(1)
PLPA.data[,c("PopNUM")][PLPA.data[,c("PopNUM")]=="PLPA_AZNE_HB"]<-as.integer(2)
PLPA.data[,c("PopNUM")][PLPA.data[,c("PopNUM")]=="PLPA_COSW_UP"]<-as.integer(3)
PLPA.data[,c("PopNUM")][PLPA.data[,c("PopNUM")]=="PLPA_UTNE_BW"]<-as.integer(4)
PLPA.data[,c("PopNUM")][PLPA.data[,c("PopNUM")]=="PLPA_UTSC_BT"]<-as.integer(5)
PLPA.data[,c("PopNUM")][PLPA.data[,c("PopNUM")]=="PLPA_UTSW_HR"]<-as.integer(6)


PLPA.RL = as.data.frame(PLPA.data[,c(55,45,49)])
PLPA.RL$PopNUM<-as.integer(PLPA.RL$PopNUM)

### 10-Days 
#PLPA.RL_10<-PLPA.RL [(PLPA.RL$Days%in%c("10")), ]
#PLPA.RL_10<-PLPA.RL_10[complete.cases(PLPA.RL_10), ]
# Compile model
#mod = stan_model("SRL_10Day.stan")
#data_list = list(alpha_mu=0,
#beta_mu=2,                  
#alpha_sigma=0,                 
#beta_sigma=.35,
#N=nrow(PLPA.RL_10),
#  n_grps = length(unique(PLPA.RL_10$PopNUM)),
#   grp = PLPA.RL_10$PopNUM,
#  y = PLPA.RL_10$S.RL)
# Sample
#PLPA_fit_SRL_10day = sampling(mod, data = data_list)

# Look at output
#print(PLPA_fit_SRL_10day, pars = c("mu", "sigma"))
#plot(PLPA_fit_SRL_10day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
#plot(PLPA_fit_SRL_10day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
#aggregate(PLPA.RL_10[,3], by=list(PoppNUM=PLPA.RL_10$PopNUM), FUN=mean)
#aggregate(PLPA.RL_10[,3], by=list(PoppNUM=PLPA.RL_10$PopNUM), FUN=sd)

### 24-Days 
PLPA.RL_24<-PLPA.RL [(PLPA.RL$Days%in%c("24")), ]
PLPA.RL_24<-PLPA.RL_24[complete.cases(PLPA.RL_24), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.28,
                 N=nrow(PLPA.RL_24),
                 n_grps = length(unique(PLPA.RL_24$PopNUM)),
                 grp = PLPA.RL_24$PopNUM,
                 y = PLPA.RL_24$S.RL)
# Sample
PLPA_fit_SRL_24day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_SRL_24day, pars = c("mu", "sigma"))
PLPA_24day_SRL_mu_plot<-plot(PLPA_fit_SRL_24day, pars = "mu",fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_24day_SRL_sigma_plot<-plot(PLPA_fit_SRL_24day, pars = "sigma",fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.RL_24[,3], by=list(PoppNUM=PLPA.RL_24$PopNUM), FUN=mean)
aggregate(PLPA.RL_24[,3], by=list(PoppNUM=PLPA.RL_24$PopNUM), FUN=sd)


### 42-Days 
PLPA.RL_42<-PLPA.RL [(PLPA.RL$Days%in%c("42")), ]
PLPA.RL_42<-PLPA.RL_42[complete.cases(PLPA.RL_42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=1.8,
                 alpha_sigma=0,
                 beta_sigma=.4,
                 N=nrow(PLPA.RL_42),
                 n_grps = length(unique(PLPA.RL_42$PopNUM)),
                 grp = PLPA.RL_42$PopNUM,
                 y = PLPA.RL_42$S.RL)
# Sample
PLPA_fit_SRL_42day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_SRL_42day, pars = c("mu", "sigma"))
PLPA_42day_SRL_mu_plot<-plot(PLPA_fit_SRL_42day, pars="mu",fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_42day_SRL_sigma_plot<-plot(PLPA_fit_SRL_42day, pars = "sigma",fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.RL_42[,3], by=list(PoppNUM=PLPA.RL_42$PopNUM), FUN=mean)
aggregate(PLPA.RL_42[,3], by=list(PoppNUM=PLPA.RL_42$PopNUM), FUN=sd)


### TOTAL ROOT LENGTH
PLPA.SumOfLength.cm. = as.data.frame(PLPA.data[,c(55,3,45)])
PLPA.SumOfLength.cm.$PopNUM<-as.integer(PLPA.SumOfLength.cm.$PopNUM)

### 10-Days 
PLPA.SumOfLength.cm._10<-PLPA.SumOfLength.cm. [(PLPA.SumOfLength.cm.$Days%in%c("10")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=3,                  
                 alpha_sigma=0,                  
                 beta_sigma=.3,
                 N=nrow(PLPA.SumOfLength.cm._10),
                 n_grps = length(unique(PLPA.SumOfLength.cm._10$PopNUM)),
                 grp = PLPA.SumOfLength.cm._10$PopNUM,
                 y = PLPA.SumOfLength.cm._10$SumOfLength.cm.)
# Sample
PLPA_fit_TRL_10day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_TRL_10day, pars = c("mu", "sigma"))
PLPA_10day_TRL_mu_plot<-plot(PLPA_fit_TRL_10day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_10day_TRL_sigma_plot<-plot(PLPA_fit_TRL_10day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))


# Check against empirical values
aggregate(PLPA.SumOfLength.cm._10[,2], by=list(PoppNUM=PLPA.SumOfLength.cm._10$PopNUM), FUN=mean)
aggregate(PLPA.SumOfLength.cm._10[,2], by=list(PoppNUM=PLPA.SumOfLength.cm._10$PopNUM), FUN=sd)

### 24-Days 
PLPA.SumOfLength.cm._24<-PLPA.SumOfLength.cm. [(PLPA.SumOfLength.cm.$Days%in%c("24")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0, 
                 beta_sigma=.3,
                 N=nrow(PLPA.SumOfLength.cm._24),
                 n_grps = length(unique(PLPA.SumOfLength.cm._24$PopNUM)),
                 grp = PLPA.SumOfLength.cm._24$PopNUM,
                 y = PLPA.SumOfLength.cm._24$SumOfLength.cm.)
# Sample
PLPA_fit_TRL_24day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_TRL_24day, pars = c("mu", "sigma"))
PLPA_24day_TRL_mu_plot<-plot(PLPA_fit_TRL_24day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_24day_TRL_sigma_plot<-plot(PLPA_fit_TRL_24day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.SumOfLength.cm._24[,2], by=list(PoppNUM=PLPA.SumOfLength.cm._24$PopNUM), FUN=mean)
aggregate(PLPA.SumOfLength.cm._24[,2], by=list(PoppNUM=PLPA.SumOfLength.cm._24$PopNUM), FUN=sd)

### 42-Days 
PLPA.SumOfLength.cm._42<-PLPA.SumOfLength.cm. [(PLPA.SumOfLength.cm.$Days%in%c("42")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0, 
                 beta_mu=7,                
                 alpha_sigma=0,               
                 beta_sigma=.3,
                 N=nrow(PLPA.SumOfLength.cm._42),
                 n_grps = length(unique(PLPA.SumOfLength.cm._42$PopNUM)),
                 grp = PLPA.SumOfLength.cm._42$PopNUM,
                 y = PLPA.SumOfLength.cm._42$SumOfLength.cm.)
# Sample
PLPA_fit_TRL_42day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_TRL_42day, pars = c("mu", "sigma"))
PLPA_42day_TRL_mu_plot<-plot(PLPA_fit_TRL_42day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_42day_TRL_sigma_plot<-plot(PLPA_fit_TRL_42day,outer_level=0.9, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.SumOfLength.cm._42[,2], by=list(PoppNUM=PLPA.SumOfLength.cm._42$PopNUM), FUN=mean)
aggregate(PLPA.SumOfLength.cm._42[,2], by=list(PoppNUM=PLPA.SumOfLength.cm._42$PopNUM), FUN=sd)

### SPECIFIC LEAF AREA 
PLPA.SLA. = as.data.frame(PLPA.data[,c(55,48,45)])
PLPA.SLA.$PopNUM<-as.integer(PLPA.SLA.$PopNUM)

### 10-Days 
#PLPA.SLA._10<-PLPA.SLA. [(PLPA.SLA.$Days%in%c("10")), ]
#PLPA.SLA._10<-PLPA.SLA._10[complete.cases(PLPA.SLA._10), ]

# Compile model
#mod = stan_model("Traits_1.stan")
#data_list = list(alpha_mu=0,
#  beta_mu=1.3,
#  alpha_sigma=0,
#  beta_sigma=.3,
#  N=nrow(PLPA.SLA._10),
#  n_grps = length(unique(PLPA.SLA._10$PopNUM)),
#  grp = PLPA.SLA._10$PopNUM,
#  y = PLPA.SLA._10$SLA.TOT)

# Sample
#PLPA_fit_SLA_10day = sampling(mod, data = data_list)

# Look at output
#print(PLPA_fit_SLA_10day, pars = c("mu", "sigma"))
#PLPA_10day_SLA_mu_plot<-plot(PLPA_fit_SLA_10day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
#PLPA_10day_SLA_sigma_plot<-plot(PLPA_fit_SLA_10day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
#aggregate(PLPA.SLA._10[,2], by=list(PoppNUM=PLPA.SLA._10$PopNUM), FUN=mean)
#aggregate(PLPA.SLA._10[,2], by=list(PoppNUM=PLPA.SLA._10$PopNUM), FUN=sd)


### 24-Days 
PLPA.SLA._24<-PLPA.SLA. [(PLPA.SLA.$Days%in%c("24")), ]
PLPA.SLA._24<-PLPA.SLA._24[complete.cases(PLPA.SLA._24), ]
#Remove row 271 (8)  because it seems incorrect - make sure to come back and check this. 
PLPA.SLA._24<-PLPA.SLA._24[-c(8),]
### Remove popNUM = 2 because it only have one value 
#PLPA.SLA._24<-PLPA.SLA._24[!(PLPA.SLA._24$PopNUM%in%c("2")),]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=6,     
                 alpha_sigma=0,     
                 beta_sigma=.2,
                 N=nrow(PLPA.SLA._24),
                 n_grps = length(unique(PLPA.SLA._24$PopNUM)),
                 grp = PLPA.SLA._24$PopNUM,
                 y = PLPA.SLA._24$SLA.TOT)
# Sample
PLPA_fit_SLA_24day = sampling(mod, data = data_list)


# Look at output
print(PLPA_fit_SLA_24day, pars = c("mu", "sigma"))
PLPA_24day_SLA_mu_plot<-plot(PLPA_fit_SLA_24day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_24day_SLA_sigma_plot<-plot(PLPA_fit_SLA_24day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.SLA._24[,2], by=list(PoppNUM=PLPA.SLA._24$PopNUM), FUN=mean)
aggregate(PLPA.SLA._24[,2], by=list(PoppNUM=PLPA.SLA._24$PopNUM), FUN=sd)

### 42-Days 
PLPA.SLA._42<-PLPA.SLA. [(PLPA.SLA.$Days%in%c("42")), ]
PLPA.SLA._42<-PLPA.SLA._42[complete.cases(PLPA.SLA._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.18,
                 N=nrow(PLPA.SLA._42),
                 n_grps = length(unique(PLPA.SLA._42$PopNUM)),
                 grp = PLPA.SLA._42$PopNUM,
                 y = PLPA.SLA._42$SLA.TOT)
# Sample
PLPA_fit_SLA_42day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_SLA_42day, pars = c("mu", "sigma"))
PLPA_42day_SLA_mu_plot<-plot(PLPA_fit_SLA_42day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_42day_SLA_sigma_plot<-plot(PLPA_fit_SLA_42day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.SLA._42[,2], by=list(PoppNUM=PLPA.SLA._42$PopNUM), FUN=mean)
aggregate(PLPA.SLA._42[,2], by=list(PoppNUM=PLPA.SLA._42$PopNUM), FUN=sd)

### HEIGHT
PLPA.HT. = as.data.frame(PLPA.data[,c(55,19,45)])
PLPA.HT.$PopNUM<-as.integer(PLPA.HT.$PopNUM)

### 10-Days 
PLPA.HT._10<-PLPA.HT. [(PLPA.HT.$Days%in%c("10")), ]
PLPA.HT._10<-PLPA.HT._10[complete.cases(PLPA.HT._10), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.2,
                 N=nrow(PLPA.HT._10),
                 n_grps = length(unique(PLPA.HT._10$PopNUM)),
                 grp = PLPA.HT._10$PopNUM,
                 y = PLPA.HT._10$HT)

# Sample
PLPA_fit_HT_10day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_HT_10day, pars = c("mu", "sigma"))
PLPA_10day_HT_mu_plot<-plot(PLPA_fit_HT_10day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_10day_HT_sigma_plot<-plot(PLPA_fit_HT_10day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.HT._10[,2], by=list(PoppNUM=PLPA.HT._10$PopNUM), FUN=mean)
aggregate(PLPA.HT._10[,2], by=list(PoppNUM=PLPA.HT._10$PopNUM), FUN=sd)

### 24-Days 
PLPA.HT._24<-PLPA.HT. [(PLPA.HT.$Days%in%c("24")), ]
PLPA.HT._24<-PLPA.HT._24[complete.cases(PLPA.HT._24), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.5,
                 alpha_sigma=0,
                 beta_sigma=.15,
                 N=nrow(PLPA.HT._24),
                 n_grps = length(unique(PLPA.HT._24$PopNUM)),
                 grp = PLPA.HT._24$PopNUM,
                 y = PLPA.HT._24$HT)
# Sample
PLPA_fit_HT_24day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_HT_24day, pars = c("mu", "sigma"))
PLPA_24day_HT_mu_plot<-plot(PLPA_fit_HT_24day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_24day_HT_sigma_plot<-plot(PLPA_fit_HT_24day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.HT._24[,2], by=list(PoppNUM=PLPA.HT._24$PopNUM), FUN=mean)
aggregate(PLPA.HT._24[,2], by=list(PoppNUM=PLPA.HT._24$PopNUM), FUN=sd)

### 42-Days 
PLPA.HT._42<-PLPA.HT. [(PLPA.HT.$Days%in%c("42")), ]
PLPA.HT._42<-PLPA.HT._42[complete.cases(PLPA.HT._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=3,  
                 alpha_sigma=0,
                 beta_sigma=.1,
                 N=nrow(PLPA.HT._42),
                 n_grps = length(unique(PLPA.HT._42$PopNUM)),
                 grp = PLPA.HT._42$PopNUM,
                 y = PLPA.HT._42$HT)
# Sample
PLPA_fit_HT_42day = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_HT_42day, pars = c("mu", "sigma"))
PLPA_42day_HT_mu_plot<-plot(PLPA_fit_HT_42day, pars = "mu", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))
PLPA_42day_HT_sigma_plot<-plot(PLPA_fit_HT_42day, pars = "sigma", fill_color=c("deepskyblue2","goldenrod3", "deepskyblue4", "olivedrab2","olivedrab4","blue2"))

# Check against empirical values
aggregate(PLPA.HT._42[,2], by=list(PoppNUM=PLPA.HT._42$PopNUM), FUN=mean)
aggregate(PLPA.HT._42[,2], by=list(PoppNUM=PLPA.HT._42$PopNUM), FUN=sd)




```

##### P. patagoinca - Total root length (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(PLPA_10day_TRL_mu_plot+xlim(0,50), PLPA_10day_TRL_sigma_plot+xlim(0,50),
             PLPA_24day_TRL_mu_plot+xlim(0,150), PLPA_24day_TRL_sigma_plot+xlim(0,150),
             PLPA_42day_TRL_mu_plot+xlim(0,600), PLPA_42day_TRL_sigma_plot+xlim(0,600),
             nrow=3 )

```

##### P. patagonica - Specific root length (cm/g)
```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
blankplot<-ggplot()
grid.arrange(blankplot,blankplot,PLPA_24day_SRL_mu_plot+xlim(0,150), PLPA_24day_SRL_sigma_plot+xlim(0,150),
             PLPA_42day_SRL_mu_plot+xlim(0,150), PLPA_42day_SRL_sigma_plot+xlim(0,150),
             nrow=3 )

```

##### P. patagonica - Specific leaf area (cm^2/g)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
blankplot<-ggplot()
grid.arrange(blankplot,blankplot,
             PLPA_24day_SLA_mu_plot+xlim(0,60), PLPA_24day_SLA_sigma_plot+xlim(0,60),
             PLPA_42day_SLA_mu_plot+xlim(0,60), PLPA_42day_SLA_sigma_plot+xlim(0,60),
             nrow=3 )

```

##### P. patagonica - Height (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(PLPA_10day_HT_mu_plot+xlim(0,4), PLPA_10day_HT_sigma_plot+xlim(0,4),
             PLPA_24day_HT_mu_plot+xlim(0,4), PLPA_24day_HT_sigma_plot+xlim(0,4),
             PLPA_42day_HT_mu_plot+xlim(0,4), PLPA_42day_HT_sigma_plot+xlim(0,4),
             nrow=3 )

```

##### P. patagonica - "Regional" SLA at 24 and 42 days 
```{r, include=FALSE, message=FALSE, warning=FALSE, error=FALSE, echo=FALSE ,results='hide'}
### SPECIFIC LEAF AREA 
PLPA.SLA. = as.data.frame(PLPA.data[,c(55,48,45)])
PLPA.SLA.$PopNUM<-as.integer(PLPA.SLA.$PopNUM)

### 24-Days 
PLPA.SLA._24<-PLPA.SLA. [(PLPA.SLA.$Days%in%c("24")), ]
PLPA.SLA._24<-PLPA.SLA._24[complete.cases(PLPA.SLA._24), ]
### Add a PopALL Column = 1 to change groups to PopALL
PLPA.SLA._24$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(PLPA.SLA._24),
                 n_grps = length(unique(PLPA.SLA._24$PopALL)),
                 grp = PLPA.SLA._24$PopALL,
                 y = PLPA.SLA._24$SLA.TOT)

# Sample
PLPA_fit_SLA_24day_PopALL = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_SLA_24day_PopALL, pars = c("mu", "sigma"))
PLPA_24day_SLA_mu_plot_PopALL<-plot(PLPA_fit_SLA_24day_PopALL, pars = "mu", fill_color="black")
PLPA_24day_SLA_sigma_plot_PopALL<-plot(PLPA_fit_SLA_24day_PopALL, pars = "sigma", fill_color=c("black"))
PLPA_24day_SLA_mu_PopALL<-plot(PLPA_fit_SLA_24day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(PLPA.SLA._24[,2], by=list(PopNUM=PLPA.SLA._24$PopALL), FUN=mean)
aggregate(PLPA.SLA._24[,2], by=list(PopNUM=PLPA.SLA._24$PopALL), FUN=sd)



### 42-Days 
PLPA.SLA._42<-PLPA.SLA. [(PLPA.SLA.$Days%in%c("42")), ]

PLPA.SLA._42<-PLPA.SLA._42[complete.cases(PLPA.SLA._42), ]
### Add a PopALL Column = 1 to change groups to PopALL
PLPA.SLA._42$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(PLPA.SLA._42),
                 n_grps = length(unique(PLPA.SLA._42$PopALL)),
                 grp = PLPA.SLA._42$PopALL,
                 y = PLPA.SLA._42$SLA.TOT)

# Sample
PLPA_fit_SLA_42day_PopALL = sampling(mod, data = data_list)

# Look at output
print(PLPA_fit_SLA_42day_PopALL, pars = c("mu", "sigma"))
PLPA_42day_SLA_mu_plot_PopALL<-plot(PLPA_fit_SLA_42day_PopALL, pars = "mu", fill_color="black")
PLPA_42day_SLA_sigma_plot_PopALL<-plot(PLPA_fit_SLA_42day_PopALL, pars = "sigma", fill_color=c("black"))
PLPA_42day_SLA_mu_PopALL<-plot(PLPA_fit_SLA_42day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(PLPA.SLA._42[,2], by=list(PoppNUM=PLPA.SLA._42$PopALL), FUN=mean)
aggregate(PLPA.SLA._42[,2], by=list(PoppNUM=PLPA.SLA._42$PopALL), FUN=sd)
```


```{r, echo=FALSE, fig.width=6.6, fig.height=3.8, warning=FALSE, error=FALSE, message=FALSE}

PLPA_24day_SLA_mu_PopALL+geom_vline(xintercept = 9.0, linetype="solid", size=.7, color="black")+
  xlim(0,30)+
  geom_vline(xintercept = 22.58, linetype="dashed", size=.9, color="deepskyblue2")+
  geom_vline(xintercept = 8.53, linetype="dashed", size=.7, color="goldenrod3")+
  geom_vline(xintercept = 15.05, linetype="dashed", size=.7, color="deepskyblue4")+
  geom_vline(xintercept = 8.52, linetype="dashed", size=.7, color="olivedrab2")+
  geom_vline(xintercept = 6.59, linetype="dashed", size=.7, color="olivedrab4")+ 
  geom_vline(xintercept = 5.34, linetype="dashed", size=.7, color="blue2")


PLPA_42day_SLA_mu_PopALL+geom_vline(xintercept = 19.35, linetype="solid", size=.7, color="black")+
  xlim(0,45)+
  geom_vline(xintercept = 31.75, linetype="dashed", size=.9, color="deepskyblue2")+
  geom_vline(xintercept = 16.01, linetype="dashed", size=.7, color="goldenrod3")+
  geom_vline(xintercept = 24.07, linetype="dashed", size=.7, color="deepskyblue4")+
  geom_vline(xintercept = 18.07, linetype="dashed", size=.7, color="olivedrab2")+
  geom_vline(xintercept = 18.76, linetype="dashed", size=.7, color="olivedrab4")+ 
  geom_vline(xintercept = 12.85, linetype="dashed", size=.7, color="blue2")

```


```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, results='hide'}
#
#
#
#
#
#
#
#
#
##
#
####
####
####
####
####
####
### VUOC

VUOC.data<-as.data.frame(SpeciesData [ which(SpeciesData$SPECIES=='VUOC'),  ])
### Duplicate POP_ID and change to PopNUM
VUOC.data$PopNUM<-VUOC.data$POP_ID

VUOC.data[,c("PopNUM")][VUOC.data[,c("PopNUM")]=="VUOC_AZSW_WB"]<-as.integer(1)
VUOC.data[,c("PopNUM")][VUOC.data[,c("PopNUM")]=="VUOC_UTSE_LM"]<-as.integer(2)
VUOC.data[,c("PopNUM")][VUOC.data[,c("PopNUM")]=="VUOC_UTSW_HR"]<-as.integer(3)
VUOC.data[,c("PopNUM")][VUOC.data[,c("PopNUM")]=="VUOC_UTC_GRC"]<-as.integer(4)
VUOC.data[,c("PopNUM")][VUOC.data[,c("PopNUM")]=="VUOC_UTSE_IM"]<-as.integer(5)


VUOC.RL = as.data.frame(VUOC.data[,c(55,45,49)])
VUOC.RL$PopNUM<-as.integer(VUOC.RL$PopNUM)

### 10-Days 
#VUOC.RL_10<-VUOC.RL [(VUOC.RL$Days%in%c("10")), ]
#VUOC.RL_10<-VUOC.RL_10[complete.cases(VUOC.RL_10), ]
# Compile model
#mod = stan_model("SRL_10Day.stan")
#data_list = list(alpha_mu=0,
#beta_mu=2,                  
#alpha_sigma=0,                 
#beta_sigma=.35,
#N=nrow(VUOC.RL_10),
#  n_grps = length(unique(VUOC.RL_10$PopNUM)),
#   grp = VUOC.RL_10$PopNUM,
#  y = VUOC.RL_10$S.RL)
# Sample
#VUOC_fit_SRL_10day = sampling(mod, data = data_list)

# Look at output
#print(VUOC_fit_SRL_10day, pars = c("mu", "sigma"))
#plot(VUOC_fit_SRL_10day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
#plot(VUOC_fit_SRL_10day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
#aggregate(VUOC.RL_10[,3], by=list(PoppNUM=VUOC.RL_10$PopNUM), FUN=mean)
#aggregate(VUOC.RL_10[,3], by=list(PoppNUM=VUOC.RL_10$PopNUM), FUN=sd)

### 24-Days 
VUOC.RL_24<-VUOC.RL [(VUOC.RL$Days%in%c("24")), ]
VUOC.RL_24<-VUOC.RL_24[complete.cases(VUOC.RL_24), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0,
                 beta_sigma=.15,
                 N=nrow(VUOC.RL_24),
                 n_grps = length(unique(VUOC.RL_24$PopNUM)),
                 grp = VUOC.RL_24$PopNUM,
                 y = VUOC.RL_24$S.RL)
# Sample
VUOC_fit_SRL_24day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_SRL_24day, pars = c("mu", "sigma"))
VUOC_24day_SRL_mu_plot<-plot(VUOC_fit_SRL_24day, pars = "mu",fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_24day_SRL_sigma_plot<-plot(VUOC_fit_SRL_24day, pars = "sigma",fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.RL_24[,3], by=list(PoppNUM=VUOC.RL_24$PopNUM), FUN=mean)
aggregate(VUOC.RL_24[,3], by=list(PoppNUM=VUOC.RL_24$PopNUM), FUN=sd)


### 42-Days 
VUOC.RL_42<-VUOC.RL [(VUOC.RL$Days%in%c("42")), ]
VUOC.RL_42<-VUOC.RL_42[complete.cases(VUOC.RL_42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.25,
                 N=nrow(VUOC.RL_42),
                 n_grps = length(unique(VUOC.RL_42$PopNUM)),
                 grp = VUOC.RL_42$PopNUM,
                 y = VUOC.RL_42$S.RL)
# Sample
VUOC_fit_SRL_42day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_SRL_42day, pars = c("mu", "sigma"))
VUOC_42day_SRL_mu_plot<-plot(VUOC_fit_SRL_42day, pars="mu",fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_42day_SRL_sigma_plot<-plot(VUOC_fit_SRL_42day, pars = "sigma",fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.RL_42[,3], by=list(PoppNUM=VUOC.RL_42$PopNUM), FUN=mean)
aggregate(VUOC.RL_42[,3], by=list(PoppNUM=VUOC.RL_42$PopNUM), FUN=sd)

### TOTAL ROOT LENGTH
VUOC.SumOfLength.cm. = as.data.frame(VUOC.data[,c(55,3,45)])
VUOC.SumOfLength.cm.$PopNUM<-as.integer(VUOC.SumOfLength.cm.$PopNUM)

### 10-Days 
VUOC.SumOfLength.cm._10<-VUOC.SumOfLength.cm. [(VUOC.SumOfLength.cm.$Days%in%c("10")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=3,                  
                 alpha_sigma=0,                  
                 beta_sigma=.3,
                 N=nrow(VUOC.SumOfLength.cm._10),
                 n_grps = length(unique(VUOC.SumOfLength.cm._10$PopNUM)),
                 grp = VUOC.SumOfLength.cm._10$PopNUM,
                 y = VUOC.SumOfLength.cm._10$SumOfLength.cm.)
# Sample
VUOC_fit_TRL_10day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_TRL_10day, pars = c("mu", "sigma"))
VUOC_10day_TRL_mu_plot<-plot(VUOC_fit_TRL_10day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_10day_TRL_sigma_plot<-plot(VUOC_fit_TRL_10day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))


# Check against empirical values
aggregate(VUOC.SumOfLength.cm._10[,2], by=list(PoppNUM=VUOC.SumOfLength.cm._10$PopNUM), FUN=mean)
aggregate(VUOC.SumOfLength.cm._10[,2], by=list(PoppNUM=VUOC.SumOfLength.cm._10$PopNUM), FUN=sd)

### 24-Days 
VUOC.SumOfLength.cm._24<-VUOC.SumOfLength.cm. [(VUOC.SumOfLength.cm.$Days%in%c("24")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2,
                 alpha_sigma=0, 
                 beta_sigma=.25,
                 N=nrow(VUOC.SumOfLength.cm._24),
                 n_grps = length(unique(VUOC.SumOfLength.cm._24$PopNUM)),
                 grp = VUOC.SumOfLength.cm._24$PopNUM,
                 y = VUOC.SumOfLength.cm._24$SumOfLength.cm.)
# Sample
VUOC_fit_TRL_24day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_TRL_24day, pars = c("mu", "sigma"))
VUOC_24day_TRL_mu_plot<-plot(VUOC_fit_TRL_24day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_24day_TRL_sigma_plot<-plot(VUOC_fit_TRL_24day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.SumOfLength.cm._24[,2], by=list(PoppNUM=VUOC.SumOfLength.cm._24$PopNUM), FUN=mean)
aggregate(VUOC.SumOfLength.cm._24[,2], by=list(PoppNUM=VUOC.SumOfLength.cm._24$PopNUM), FUN=sd)

### 42-Days 
VUOC.SumOfLength.cm._42<-VUOC.SumOfLength.cm. [(VUOC.SumOfLength.cm.$Days%in%c("42")), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0, 
                 beta_mu=7,                
                 alpha_sigma=0,               
                 beta_sigma=.33,
                 N=nrow(VUOC.SumOfLength.cm._42),
                 n_grps = length(unique(VUOC.SumOfLength.cm._42$PopNUM)),
                 grp = VUOC.SumOfLength.cm._42$PopNUM,
                 y = VUOC.SumOfLength.cm._42$SumOfLength.cm.)
# Sample
VUOC_fit_TRL_42day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_TRL_42day, pars = c("mu", "sigma"))
VUOC_42day_TRL_mu_plot<-plot(VUOC_fit_TRL_42day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_42day_TRL_sigma_plot<-plot(VUOC_fit_TRL_42day,outer_level=0.9, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.SumOfLength.cm._42[,2], by=list(PoppNUM=VUOC.SumOfLength.cm._42$PopNUM), FUN=mean)
aggregate(VUOC.SumOfLength.cm._42[,2], by=list(PoppNUM=VUOC.SumOfLength.cm._42$PopNUM), FUN=sd)

### SPECIFIC LEAF AREA 
VUOC.SLA. = as.data.frame(VUOC.data[,c(55,48,45)])
VUOC.SLA.$PopNUM<-as.integer(VUOC.SLA.$PopNUM)

### 10-Days 
#VUOC.SLA._10<-VUOC.SLA. [(VUOC.SLA.$Days%in%c("10")), ]
#VUOC.SLA._10<-VUOC.SLA._10[complete.cases(VUOC.SLA._10), ]

# Compile model
#mod = stan_model("Traits_1.stan")
#data_list = list(alpha_mu=0,
#  beta_mu=1.3,
#  alpha_sigma=0,
#  beta_sigma=.3,
#  N=nrow(VUOC.SLA._10),
#  n_grps = length(unique(VUOC.SLA._10$PopNUM)),
#  grp = VUOC.SLA._10$PopNUM,
#  y = VUOC.SLA._10$SLA.TOT)

# Sample
#VUOC_fit_SLA_10day = sampling(mod, data = data_list)

# Look at output
#print(VUOC_fit_SLA_10day, pars = c("mu", "sigma"))
#VUOC_10day_SLA_mu_plot<-plot(VUOC_fit_SLA_10day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", deepskyblue2"))
#VUOC_10day_SLA_sigma_plot<-plot(VUOC_fit_SLA_10day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
#aggregate(VUOC.SLA._10[,2], by=list(PoppNUM=VUOC.SLA._10$PopNUM), FUN=mean)
#aggregate(VUOC.SLA._10[,2], by=list(PoppNUM=VUOC.SLA._10$PopNUM), FUN=sd)


### 24-Days 
VUOC.SLA._24<-VUOC.SLA. [(VUOC.SLA.$Days%in%c("24")), ]
VUOC.SLA._24<-VUOC.SLA._24[complete.cases(VUOC.SLA._24), ]
#Remove row 271 (8)  because it seems incorrect - make sure to come back and check this. 
VUOC.SLA._24<-VUOC.SLA._24[-c(8),]
### Remove popNUM = 2 because it only have one value 
#VUOC.SLA._24<-VUOC.SLA._24[!(VUOC.SLA._24$PopNUM%in%c("2")),]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=6,     
                 alpha_sigma=0,     
                 beta_sigma=.22,
                 N=nrow(VUOC.SLA._24),
                 n_grps = length(unique(VUOC.SLA._24$PopNUM)),
                 grp = VUOC.SLA._24$PopNUM,
                 y = VUOC.SLA._24$SLA.TOT)
# Sample
VUOC_fit_SLA_24day = sampling(mod, data = data_list)


# Look at output
print(VUOC_fit_SLA_24day, pars = c("mu", "sigma"))
VUOC_24day_SLA_mu_plot<-plot(VUOC_fit_SLA_24day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_24day_SLA_sigma_plot<-plot(VUOC_fit_SLA_24day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.SLA._24[,2], by=list(PoppNUM=VUOC.SLA._24$PopNUM), FUN=mean)
aggregate(VUOC.SLA._24[,2], by=list(PoppNUM=VUOC.SLA._24$PopNUM), FUN=sd)

### 42-Days 
VUOC.SLA._42<-VUOC.SLA. [(VUOC.SLA.$Days%in%c("42")), ]
VUOC.SLA._42<-VUOC.SLA._42[complete.cases(VUOC.SLA._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.18,
                 N=nrow(VUOC.SLA._42),
                 n_grps = length(unique(VUOC.SLA._42$PopNUM)),
                 grp = VUOC.SLA._42$PopNUM,
                 y = VUOC.SLA._42$SLA.TOT)
# Sample
VUOC_fit_SLA_42day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_SLA_42day, pars = c("mu", "sigma"))
VUOC_42day_SLA_mu_plot<-plot(VUOC_fit_SLA_42day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_42day_SLA_sigma_plot<-plot(VUOC_fit_SLA_42day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.SLA._42[,2], by=list(PoppNUM=VUOC.SLA._42$PopNUM), FUN=mean)
aggregate(VUOC.SLA._42[,2], by=list(PoppNUM=VUOC.SLA._42$PopNUM), FUN=sd)

### HEIGHT
VUOC.HT. = as.data.frame(VUOC.data[,c(55,19,45)])
VUOC.HT.$PopNUM<-as.integer(VUOC.HT.$PopNUM)

### 10-Days 
VUOC.HT._10<-VUOC.HT. [(VUOC.HT.$Days%in%c("10")), ]
VUOC.HT._10<-VUOC.HT._10[complete.cases(VUOC.HT._10), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.2,
                 alpha_sigma=0,
                 beta_sigma=.2,
                 N=nrow(VUOC.HT._10),
                 n_grps = length(unique(VUOC.HT._10$PopNUM)),
                 grp = VUOC.HT._10$PopNUM,
                 y = VUOC.HT._10$HT)

# Sample
VUOC_fit_HT_10day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_HT_10day, pars = c("mu", "sigma"))
VUOC_10day_HT_mu_plot<-plot(VUOC_fit_HT_10day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_10day_HT_sigma_plot<-plot(VUOC_fit_HT_10day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.HT._10[,2], by=list(PoppNUM=VUOC.HT._10$PopNUM), FUN=mean)
aggregate(VUOC.HT._10[,2], by=list(PoppNUM=VUOC.HT._10$PopNUM), FUN=sd)

### 24-Days 
VUOC.HT._24<-VUOC.HT. [(VUOC.HT.$Days%in%c("24")), ]
VUOC.HT._24<-VUOC.HT._24[complete.cases(VUOC.HT._24), ]

# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=2.5,
                 alpha_sigma=0,
                 beta_sigma=.15,
                 N=nrow(VUOC.HT._24),
                 n_grps = length(unique(VUOC.HT._24$PopNUM)),
                 grp = VUOC.HT._24$PopNUM,
                 y = VUOC.HT._24$HT)
# Sample
VUOC_fit_HT_24day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_HT_24day, pars = c("mu", "sigma"))
VUOC_24day_HT_mu_plot<-plot(VUOC_fit_HT_24day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_24day_HT_sigma_plot<-plot(VUOC_fit_HT_24day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.HT._24[,2], by=list(PoppNUM=VUOC.HT._24$PopNUM), FUN=mean)
aggregate(VUOC.HT._24[,2], by=list(PoppNUM=VUOC.HT._24$PopNUM), FUN=sd)

### 42-Days 
VUOC.HT._42<-VUOC.HT. [(VUOC.HT.$Days%in%c("42")), ]
VUOC.HT._42<-VUOC.HT._42[complete.cases(VUOC.HT._42), ]
# Compile model
mod = stan_model("Traits_1.stan")
data_list = list(alpha_mu=0,
                 beta_mu=3,  
                 alpha_sigma=0,
                 beta_sigma=.1,
                 N=nrow(VUOC.HT._42),
                 n_grps = length(unique(VUOC.HT._42$PopNUM)),
                 grp = VUOC.HT._42$PopNUM,
                 y = VUOC.HT._42$HT)
# Sample
VUOC_fit_HT_42day = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_HT_42day, pars = c("mu", "sigma"))
VUOC_42day_HT_mu_plot<-plot(VUOC_fit_HT_42day, pars = "mu", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))
VUOC_42day_HT_sigma_plot<-plot(VUOC_fit_HT_42day, pars = "sigma", fill_color=c("orangered2","olivedrab2", "olivedrab4", "forestgreen", "deepskyblue2"))

# Check against empirical values
aggregate(VUOC.HT._42[,2], by=list(PoppNUM=VUOC.HT._42$PopNUM), FUN=mean)
aggregate(VUOC.HT._42[,2], by=list(PoppNUM=VUOC.HT._42$PopNUM), FUN=sd)




```

##### V. octoflora - Total root length (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(VUOC_10day_TRL_mu_plot+xlim(0,50), VUOC_10day_TRL_sigma_plot+xlim(0,50),
             VUOC_24day_TRL_mu_plot+xlim(0,170), VUOC_24day_TRL_sigma_plot+xlim(0,170),
             VUOC_42day_TRL_mu_plot+xlim(0,1200), VUOC_42day_TRL_sigma_plot+xlim(0,1200),
             nrow=3 )

```

##### V. octoflora - Specific root length (cm/g)
```{r, echo=FALSE, warning=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
blankplot<-ggplot()
grid.arrange(blankplot,blankplot,VUOC_24day_SRL_mu_plot+xlim(0,100), VUOC_24day_SRL_sigma_plot+xlim(0,100),
             VUOC_42day_SRL_mu_plot+xlim(0,100), VUOC_42day_SRL_sigma_plot+xlim(0,100),
             nrow=3 )

```

##### V. octoflora - Specific leaf area (cm^2/g)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
blankplot<-ggplot()
grid.arrange(blankplot,blankplot,
             VUOC_24day_SLA_mu_plot+xlim(0,60), VUOC_24day_SLA_sigma_plot+xlim(0,60),
             VUOC_42day_SLA_mu_plot+xlim(0,60), VUOC_42day_SLA_sigma_plot+xlim(0,60),
             nrow=3 )

```

##### V. octoflora - Height (cm)
```{r, echo=FALSE, message=FALSE, error=FALSE, warning=FALSE, fig.width=6.6, fig.height=3.8}
### Plots of mu and sigma at all three time points together:
grid.arrange(VUOC_10day_HT_mu_plot+xlim(0,4), VUOC_10day_HT_sigma_plot+xlim(0,4),
             VUOC_24day_HT_mu_plot+xlim(0,4), VUOC_24day_HT_sigma_plot+xlim(0,4),
             VUOC_42day_HT_mu_plot+xlim(0,4), VUOC_42day_HT_sigma_plot+xlim(0,4),
             nrow=3 )

```

##### V. octoflora - "Regional" SLA at 24 and 42 days 
```{r, include=FALSE, message=FALSE, warning=FALSE, error=FALSE, echo=FALSE ,results='hide'}
### SPECIFIC LEAF AREA 
VUOC.SLA. = as.data.frame(VUOC.data[,c(55,48,45)])
VUOC.SLA.$PopNUM<-as.integer(VUOC.SLA.$PopNUM)

### 24-Days 
VUOC.SLA._24<-VUOC.SLA. [(VUOC.SLA.$Days%in%c("24")), ]
VUOC.SLA._24<-VUOC.SLA._24[complete.cases(VUOC.SLA._24), ]
### Add a PopALL Column = 1 to change groups to PopALL
VUOC.SLA._24$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(VUOC.SLA._24),
                 n_grps = length(unique(VUOC.SLA._24$PopALL)),
                 grp = VUOC.SLA._24$PopALL,
                 y = VUOC.SLA._24$SLA.TOT)

# Sample
VUOC_fit_SLA_24day_PopALL = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_SLA_24day_PopALL, pars = c("mu", "sigma"))
VUOC_24day_SLA_mu_plot_PopALL<-plot(VUOC_fit_SLA_24day_PopALL, pars = "mu", fill_color="black")
VUOC_24day_SLA_sigma_plot_PopALL<-plot(VUOC_fit_SLA_24day_PopALL, pars = "sigma", fill_color=c("black"))
VUOC_24day_SLA_mu_PopALL<-plot(VUOC_fit_SLA_24day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(VUOC.SLA._24[,2], by=list(PopNUM=VUOC.SLA._24$PopALL), FUN=mean)
aggregate(VUOC.SLA._24[,2], by=list(PopNUM=VUOC.SLA._24$PopALL), FUN=sd)



### 42-Days 
VUOC.SLA._42<-VUOC.SLA. [(VUOC.SLA.$Days%in%c("42")), ]

VUOC.SLA._42<-VUOC.SLA._42[complete.cases(VUOC.SLA._42), ]
### Add a PopALL Column = 1 to change groups to PopALL
VUOC.SLA._42$PopALL<-1

# Compile model
mod = stan_model("Traits_2.stan")
data_list = list(alpha_mu=0,
                 beta_mu=.6,
                 alpha_sigma=0,
                 beta_sigma=.20,
                 N=nrow(VUOC.SLA._42),
                 n_grps = length(unique(VUOC.SLA._42$PopALL)),
                 grp = VUOC.SLA._42$PopALL,
                 y = VUOC.SLA._42$SLA.TOT)

# Sample
VUOC_fit_SLA_42day_PopALL = sampling(mod, data = data_list)

# Look at output
print(VUOC_fit_SLA_42day_PopALL, pars = c("mu", "sigma"))
VUOC_42day_SLA_mu_plot_PopALL<-plot(VUOC_fit_SLA_42day_PopALL, pars = "mu", fill_color="black")
VUOC_42day_SLA_sigma_plot_PopALL<-plot(VUOC_fit_SLA_42day_PopALL, pars = "sigma", fill_color=c("black"))


VUOC_42day_SLA_mu_PopALL<-plot(VUOC_fit_SLA_42day_PopALL, pars="mu", show_density = TRUE, ci_level = .95,outer_level=1, fill_color = "grey95")


# Check against empirical values
aggregate(VUOC.SLA._42[,2], by=list(PoppNUM=VUOC.SLA._42$PopALL), FUN=mean)
aggregate(VUOC.SLA._42[,2], by=list(PoppNUM=VUOC.SLA._42$PopALL), FUN=sd)
```


