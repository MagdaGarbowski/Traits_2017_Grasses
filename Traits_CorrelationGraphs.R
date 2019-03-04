### Analyses Code for Trait Data 2017 
### 03-03-2019 

library(lattice)
SpeciesData<-read.csv("HECO_ELTR_PLPA_VUOC.csv")
# Keep relevant columns 


# Ideally, I want to work with all of th species and all of the time points so I can compare:
# (1) within a species, among time points
# (2) within a timepoint, among populations 
# (3) within a timepoint, among species


### Can I add regression lines for each "group"? What about over all? 

### Function for correlation graphs with r and p values 

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  if(p<0.05) txt2<-paste("* ", txt2, sep = "")
  text(0.5, 0.4, txt2)
}


### All HECO 
HECO.data<-as.data.frame(SpeciesData [ which(SpeciesData$SPECIES=='HECO'),  ])
HECO.data$Days<-as.factor(HECO.data$Days)
HECO.data.complete<-na.omit(HECO.data)
HECO.data.complete2<-HECO.data.complete[c(3,4,6,7,11,12,13,14,17)]

### 10 Day values 
HECO.data.10day<-HECO.data.complete[HECO.data.complete$Days%in%c("10"),]
HECO.data.10.day.complete2<-HECO.data.10day[c(3,4,6,7,11,12,13,14,17)]
### 24 Day values 
HECO.data.24day<-HECO.data.complete[HECO.data.complete$Days%in%c("24"),]
HECO.data.24.day.complete2<-HECO.data.24day[c(3,4,6,7,11,12,13,14,17)]
### 48 Day values 
HECO.data.48day<-HECO.data.complete[HECO.data.complete$Days%in%c("42"),]
HECO.data.48.day.complete2<-HECO.data.48day[c(3,4,6,7,11,12,13,14,17)]


### DIFFERENT POPULATIONS - BUT LUMPING ALL TIME POINTS TOGETHER 

HECO.data.HECO_AZNC_KV<-HECO.data.complete[HECO.data.complete$POP_ID%in%c("HECO_AZNC_KV"),]
HECO.data.HECO_AZNC_KV.complete2<-HECO.data.HECO_AZNC_KV[c(3,4,6,7,11,12,13,14,17)]
# 
HECO.data.HECO_UTC_KRR<-HECO.data.complete[HECO.data.complete$POP_ID%in%c("HECO_UTC_KRR"),]
HECO.data.HECO_UTC_KRR.complete2<-HECO.data.HECO_UTC_KRR[c(3,4,6,7,11,12,13,14,17)]
#
HECO.data.HECO_UTC_TMR<-HECO.data.complete[HECO.data.complete$POP_ID%in%c("HECO_UTC_TMR"),]
HECO.data.HECO_UTC_TMR.complete2<-HECO.data.HECO_UTC_TMR[c(3,4,6,7,11,12,13,14,17)]
#
HECO.data.HECO_UTEC_GR<-HECO.data.complete[HECO.data.complete$POP_ID%in%c("HECO_UTEC_GR"),]
HECO.data.HECO_UTEC_GR.complete2<-HECO.data.24day[c(3,4,6,7,11,12,13,14,17)]
# 
HECO.data.HECO_UTSW_SP<-HECO.data.complete[HECO.data.complete$POP_ID%in%c("HECO_UTSW_SP"),]
HECO.data.HECO_UTSW_SP.complete2<-HECO.data.HECO_UTSW_SP[c(3,4,6,7,11,12,13,14,17)]

#### Time point correlations
pairs(HECO.data.HECO_AZNC_KV.complete2, upper.panel = panel.cor, main="HECO_AZNC_KV")

#### Population correlations 

pairs(HECO.data.10.day.complete2, upper.panel = panel.cor, main="HECO_10Day")
pairs(HECO.data.24.day.complete2, upper.panel = panel.cor, main="HECO_24Day")
pairs(HECO.data.42.day.complete2, upper.panel = panel.cor, main="HECO_42Day")

### All ELTR 
ELTR.data<-as.data.frame(SpeciesData [ which(SpeciesData$SPECIES=='ELTR'),  ])
ELTR.data$Days<-as.factor(ELTR.data$Days)
ELTR.data.complete<-na.omit(ELTR.data)
ELTR.data.complete2<-ELTR.data.complete[c(3,4,6,7,11,12,13,14,17)]

### 10 Day values 
ELTR.data.10day<-ELTR.data.complete[ELTR.data.complete$Days%in%c("10"),]
ELTR.data.10.day.complete2<-ELTR.data.10day[c(3,4,6,7,11,12,13,14,17)]
### 24 Day values 
ELTR.data.24day<-ELTR.data.complete[ELTR.data.complete$Days%in%c("24"),]
ELTR.data.24.day.complete2<-ELTR.data.24day[c(3,4,6,7,11,12,13,14,17)]
### 48 Day values 
ELTR.data.42day<-ELTR.data.complete[ELTR.data.complete$Days%in%c("42"),]
ELTR.data.42.day.complete2<-ELTR.data.42day[c(3,4,6,7,11,12,13,14,17)]


# ALL POPULATIONS BUT LUMPING TIME POINTS TOGETHER
ELTR.data.ELTR_NMNW_C_<-ELTR.data.complete[ELTR.data.complete$POP_ID%in%c("ELTR_NMNW_C_"),]
ELTR.data.ELTR_NMNW_C_.complete2<-ELTR.data.ELTR_NMNW_C_[c(3,4,6,7,11,12,13,14,17)]
### 
ELTR.data.ELTR_UTC_TMC<-ELTR.data.complete[ELTR.data.complete$POP_ID%in%c("ELTR_UTC_TMC"),]
ELTR.data.ELTR_UTC_TMC.complete2<-ELTR.data.ELTR_UTC_TMC[c(3,4,6,7,11,12,13,14,17)]
### 
ELTR.data.ELTR_UTSW_DN<-ELTR.data.complete[ELTR.data.complete$POP_ID%in%c("ELTR_UTSW_DN"),]
ELTR.data.ELTR_UTSW_DN.complete2<-ELTR.data.ELTR_UTSW_DN[c(3,4,6,7,11,12,13,14,17)]
### 
ELTR.data.ELTR_UTSW_SPR<-ELTR.data.complete[ELTR.data.complete$POP_ID%in%c("ELTR_UTSW_SP"),]
ELTR.data.ELTR_UTSW_SP.complete2<-ELTR.data.24day[c(3,4,6,7,11,12,13,14,17)]


#### Time point correlations
pairs(ELTR.data.complete2, upper.panel = panel.cor, main = "ALL ELTR")
pairs(ELTR.data.10.day.complete2, upper.panel = panel.cor, main="ELTR_10Day")
pairs(ELTR.data.24.day.complete2, upper.panel = panel.cor, main="ELTR_24Day")
pairs(ELTR.data.42.day.complete2, upper.panel = panel.cor, main="ELTR_42Day")

### Population correlations 

pairs(ELTR.data.ELTR_NMNW_C_.complete2, upper.panel = panel.cor, main="ELTR_NMNW_C_")
pairs(ELTR.data.ELTR_UTC_TMC.complete2, upper.panel = panel.cor, main="ELTR_UTC_TMC")
pairs(ELTR.data.ELTR_UTSW_DN.complete2, upper.panel = panel.cor, main="ELTR_UTSW_DN")
pairs(ELTR.data.ELTR_UTSW_SP.complete2, upper.panel = panel.cor, main="ELTR_UTSW_SP")


### Complete correlations for VUOC 
### Rename as correlations code or something...


