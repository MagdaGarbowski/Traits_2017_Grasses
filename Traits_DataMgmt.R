### Data Mgmt Code for Trait Data 2017 
### 03-03-2019 
### Largely working from: HECO_HEVI_PLPA_VUOC_Bayes_PointEstimateGraphs.Rmd

#### H. comata point estimates - mean and variance - 3/21/2018

library(rstan)
library(plyr)
library(data.table)
library(ggplot2)
library(gridExtra)

### Load trait data for HECO, HEVI, PLPA, VUOC 

HECODATA<-read.csv("/Users/MagdaGarbowski/Traits/Data_Query_HECO.csv")
ELTRDATA<-read.csv("/Users/MagdaGarbowski/Traits/Data_Query_ELTR.csv")
VUOCDATA<-read.csv("/Users/MagdaGarbowski/Traits/Data_Query_VUOC.csv")
PLPADATA<-read.csv("/Users/MagdaGarbowski/Traits/Data_Query_PLPA.csv")
### Remove "late" measurements from PLPA 
PLPADATA<-PLPADATA [(PLPADATA$NOTES%in%c("Early")), ]       #Select out just early species 

SpeciesData<-rbind(HECODATA, ELTRDATA, VUOCDATA, PLPADATA)
head(SpeciesData)

### Create POP_ID Column
SpeciesData$POP_ID <- substr(SpeciesData$SAMPLE_ID,1,12) # Create population ID column
SpeciesData$H_num<-as.integer(SpeciesData$H_num)
### Create a column for days so that x-axis on graph says days not harvests 
SpeciesData$Days<-as.factor(SpeciesData$H_num) # Duplicate H_num column
levels(SpeciesData$Days)<-c("10","24","42", "84") # Replace with "Day" values

### Create Species Column 
SpeciesData$SPECIES <- substr(SpeciesData$POP_ID,1,4)

### Create columns for SLA, SRL, Total Root Length, RMR, Height, Total Weight, LDMC

###SLA TOT
SpeciesData$TOTALEAFWEIGHT<- rowSums(SpeciesData[,c("LWD", "LWD_A")], na.rm=TRUE)
SpeciesData$SLA.TOT<-(SpeciesData$TOTALEAFWEIGHT/SpeciesData$Total.Of.PROJ_AREA_AG)*10000 # ALL TOGETHER
str(SpeciesData$SLA.TOT)
is.na(SpeciesData$SLA.TOT) <- !SpeciesData$SLA.TOT

### SRL 
### Create a column for S.RL = SumofSumOfLength.cm./RWD
SpeciesData$S.RL<-(SpeciesData$RWD/SpeciesData$SumOfLength.cm.)*1000000
str(SpeciesData$S.RL)
is.na(SpeciesData$S.RL) <- !SpeciesData$S.RL

### RMR
SpeciesData$TotWD<-rowSums(SpeciesData[,c("RWD","LWD","LWD_A","SWD","CWD")],na.rm=TRUE) 
SpeciesData$RMR<-SpeciesData$RWD/SpeciesData$TotWD
SpeciesData[,c("RMR")][SpeciesData[,c("RMR")]==0]<-NA   ### Make 0s and 1s NA

### LDMC
SpeciesData$TotLDW<-rowSums(SpeciesData[,c("LWD","LWD_A","SWD")],na.rm=TRUE) 
SpeciesData$TotLWS<-rowSums(SpeciesData[,c("LWS","LWS_A","SWS")],na.rm=TRUE) 
SpeciesData$LDMC<-SpeciesData$TotLDW / SpeciesData$TotLWS
str(SpeciesData$LDMC)
is.na(SpeciesData$LDMC) <- !SpeciesData$LDMC

SpeciesData[,c("LDMC")][SpeciesData[,c("LDMC")]==0]<-NA


### Keep relevant columns 
SpeciesData<- as.data.frame(SpeciesData[, c(1,3,6,18,19,31,44,45,46,48:54)])

### Write complete dataset to CSV 
write.csv(SpeciesData, "HECO_ELTR_PLPA_VUOC.csv")
