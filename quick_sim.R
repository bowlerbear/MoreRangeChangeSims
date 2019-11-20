set.seed(3)

#trend
type=21
getParams(parameterDF,mysimNu=type)
Obs <- replicate(300,getSims(),simplify=FALSE)
save(Obs,file=paste0("output/Obs_",type,"_",systemTime,".RData"))

myModel <- "basic_odm_fixedYear_LL.txt"
out<- lapply(Obs,function(x)getModels(x,bias="standard",myModel))
save(out,file=paste0("output/out_standard",type,"_",systemTime,".RData"))
outRaw <- ldply(Obs,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_standard",type,"_",systemTime,".RData"))

type=23
myModel <- "basic_odm_randomwalk_LL.txt"
out<- lapply(Obs,function(x)getModels(x,bias="standard",myModel))
save(out,file=paste0("output/out_standard",type,"_",systemTime,".RData"))
outRaw <- ldply(Obs,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_standard",type,"_",systemTime,".RData"))

type=25
myModel <- "basic_odm_dynamic_LL.txt"
out<- lapply(Obs,function(x)getModels(x,bias="standard",myModel))
save(out,file=paste0("output/out_standard",type,"_",systemTime,".RData"))
outRaw <- ldply(Obs,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_standard",type,"_",systemTime,".RData"))

#autocorrelation scenarios
type=21
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_autocorrelation.R')

myModel <- "basic_odm_fixedYear_LL.txt"
(1)
runSims <- function(Obs,propSites=1,autoProb=0.5){
  #reduce visits in non-atlas years
  Obs_auto <-applyAutocorrelation(Obs,propSites,autoProb)
  return(Obs_auto)
}
Obs_auto <- lapply(Obs,runSims)
#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_Autocorrelation100",type,"_",systemTime,".RData"))
#fit standard OD model
out <- lapply(Obs_auto,function(x)getModels(x,bias="Autocorrelation100",myModel))
save(out,file=paste0("output/out_Autocorrelation100",type,"_",systemTime,".RData"))

(2)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_autocorrelation.R')

runSims <- function(Obs,propSites=1){
  #reduce visits in non-atlas years
  Obs_auto <-applyFullAutocorrelation(Obs,propSites)
  return(Obs_auto)
}

Obs_auto <- lapply(Obs,runSims)
#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_FullAutocorrelation100",type,"_",systemTime,".RData"))
#fit standard OD model
out <- lapply(Obs_auto,function(x)getModels(x,bias="FullAutocorrelation100",myModel))
save(out,file=paste0("output/out_FullAutocorrelation100",type,"_",systemTime,".RData"))

(3)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_pulsedvisits.R')

runSims <- function(Obs,nStandardVisits=3){
  #reduce visits in non-atlas years
  Obs_pulsed <- pulsedVisits(Obs, nStandardVisits)
  #return summary data
  return(Obs_pulsed)
}

Obs_auto <- lapply(Obs,runSims)
#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_Pulsed",type,"_",systemTime,".RData"))
#fit OD models
out <- lapply(Obs_auto,function(x)getModels(x,bias="Pulsed",myModel))
save(out,file=paste0("output/out_Pulsed",type,"_",systemTime,".RData"))

(4)  
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_atlassites.R')

runSims <- function(Obs,nStandardSites=25){
  #reduce visits in non-atlas years
  Obs_spread <- reduceSites(Obs,nStandardSites)
  #return summary data
  return(Obs_spread)
}
Obs_auto <- lapply(Obs,runSims)

#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_Spread",type,"_",systemTime,".RData"))
#fit OD models
out <- lapply(Obs_auto,function(x)getModels(x,bias="Spread",myModel))
save(out,file=paste0("output/out_Spread",type,"_",systemTime,".RData"))

(5)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_pulsedvisits.R')

runSims <- function(Obs,lowProb=0.25){
  #reduce visits in non-atlas years
  Obs_spread <- pulsedSpreadVisits(Obs,SCovariate,lowProb)
  #return summary data
  return(Obs_spread)
}

Obs_auto <- lapply(Obs,runSims)

#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_SpreadLower",type,"_",systemTime,".RData"))
#fit OD models
out_spread <- lapply(Obs_auto,function(x)getModels(x,bias="SpreadLower",myModel))
save(out,file=paste0("output/out_SpreadLower",type,"_",systemTime,".RData"))

(6)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_frequencyvisits.R')

runSims <- function(Obs,type="static"){
  Obs_reduced <-applyReducedVisits(Obs,type)
  return(Obs_reduced)
}

Obs_auto <- lapply(Obs,runSims)

#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_ReducedQ",type,"_",systemTime,".RData"))
#fit standard OD model
out <- lapply(Obs_auto,function(x)getModels(x,bias="ReducedQ",myModel))
save(out,file=paste0("output/out_ReducedQ",type,"_",systemTime,".RData"))


##############################################################################################

#null trend
type=22
getParams(parameterDF,mysimNu=type)
Obs <- replicate(300,getSims(),simplify=FALSE)
save(Obs,file=paste0("output/Obs_",type,"_",systemTime,".RData"))

myModel <- "basic_odm_fixedYear_LL.txt"
out<- lapply(Obs,function(x)getModels(x,bias="standard",myModel))
save(out,file=paste0("output/out_standard",type,"_",systemTime,".RData"))
outRaw <- ldply(Obs,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_standard",type,"_",systemTime,".RData"))


type=24
myModel <- "basic_odm_randomwalk_LL.txt"
out<- lapply(Obs,function(x)getModels(x,bias="standard",myModel))
save(out,file=paste0("output/out_standard",type,"_",systemTime,".RData"))
outRaw <- ldply(Obs,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_standard",type,"_",systemTime,".RData"))

type=26
myModel <- "basic_odm_dynamic_LL.txt"
out<- lapply(Obs,function(x)getModels(x,bias="standard",myModel))
save(out,file=paste0("output/out_standard",type,"_",systemTime,".RData"))
outRaw <- ldply(Obs,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_standard",type,"_",systemTime,".RData"))

#autocorrelation scenarios
type=22
myModel <- "basic_odm_fixedYear_LL.txt"
(1)
runSims <- function(Obs,propSites=0.5,autoProb=0.5){
  #reduce visits in non-atlas years
  Obs_auto <-applyAutocorrelation(Obs,propSites,autoProb)
  return(Obs_auto)
}
Obs_auto <- lapply(Obs,runSims)
#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_Autocorrelation",type,"_",systemTime,".RData"))
#fit standard OD model
out <- lapply(Obs_auto,function(x)getModels(x,bias="Autocorrelation",myModel))
save(out,file=paste0("output/out_Autocorrelation",type,"_",systemTime,".RData"))
#fix
outFixed <- lapply(Obs_auto,function(x)fixAutocorrelation(x))
save(outFixed,file=paste0("output/outFixed_Autocorrelation",type,"_",systemTime,".RData"))



(2)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_autocorrelation.R')

runSims <- function(Obs,propSites=0.5){
  #reduce visits in non-atlas years
  Obs_auto <-applyFullAutocorrelation(Obs,propSites)
  return(Obs_auto)
}

Obs_auto <- lapply(Obs,runSims)
#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_FullAutocorrelation",type,"_",systemTime,".RData"))
#fit standard OD model
out <- lapply(Obs_auto,function(x)getModels(x,bias="FullAutocorrelation",myModel))
save(out,file=paste0("output/out_FullAutocorrelation",type,"_",systemTime,".RData"))
outFixed <- lapply(Obs_auto,function(x)fixFullAutocorrelation(x))
save(outFixed,file=paste0("output/outFixed_FullAutocorrelation",type,"_",systemTime,".RData"))


(3)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_pulsedvisits.R')

runSims <- function(Obs,nStandardVisits=3){
  #reduce visits in non-atlas years
  Obs_pulsed <- pulsedVisits(Obs, nStandardVisits)
  #return summary data
  return(Obs_pulsed)
}

Obs_auto <- lapply(Obs,runSims)
#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_Pulsed",type,"_",systemTime,".RData"))
#fit OD models
out <- lapply(Obs_auto,function(x)getModels(x,bias="Pulsed",myModel))
save(out,file=paste0("output/out_Pulsed",type,"_",systemTime,".RData"))

(4) 
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_atlassites.R')

runSims <- function(Obs,nStandardSites=25){
  #reduce visits in non-atlas years
  Obs_spread <- reduceSites(Obs,nStandardSites)
  #return summary data
  return(Obs_spread)
}
Obs_auto <- lapply(Obs,runSims)

#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_Spread",type,"_",systemTime,".RData"))
#fit OD models
out <- lapply(Obs_auto,function(x)getModels(x,bias="Spread",myModel))
save(out,file=paste0("output/out_Spread",type,"_",systemTime,".RData"))

(5)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_pulsedvisits.R')

runSims <- function(Obs,lowProb=0.25){
  #reduce visits in non-atlas years
  Obs_spread <- pulsedSpreadVisits(Obs,SCovariate,lowProb)
  #return summary data
  return(Obs_spread)
}

Obs_auto <- lapply(Obs,runSims)

#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_SpreadLower",type,"_",systemTime,".RData"))
#fit OD models
out_spread <- lapply(Obs_auto,function(x)getModels(x,bias="SpreadLower",myModel))
save(out,file=paste0("output/out_SpreadLower",type,"_",systemTime,".RData"))







(6)
source('C:/Users/db40fysa/Nextcloud/sMon-Analyses/MoreRangeChangeSims/R/scenario_frequencyvisits.R')

runSims <- function(Obs,type="static"){
  Obs_reduced <-applyReducedVisits(Obs,type)
  return(Obs_reduced)
}

Obs_auto <- lapply(Obs,runSims)

#raw data analysis
outRaw<- ldply(Obs_auto,function(x)rawdataAnalysis(x))
save(outRaw,file=paste0("output/outRaw_ReducedQ",type,"_",systemTime,".RData"))

#fit standard OD model
out <- lapply(Obs_auto,function(x)getModels(x,bias="ReducedQ",myModel))
save(out,file=paste0("output/out_ReducedQ",type,"_",systemTime,".RData"))


#####################################
#type 21



######################################

#options (1)


# Likelihood model   
for(i in 1:R){            
  z[i] ~ dbern(psi)
  p.occ[i] <- z[i]*p
  counts[i] ~ dbinom(p.occ[i],camdays[i])
}

#Priors
psi ~ dunif(0,1)
p ~ dunif(0,1)

# Addition to the likelihood model
logit(p[i]) <- alpha + betaM * Margin[i]

# Additional priors
alpha ~ dunif(-20,20)
betaM ~ dunif(-20,20)

# Addition to the likelihood model
p[i] ~ dbeta(prior.a,prior.b)

# Additional priors
prior.avg ~ dbeta(1,1)
prior.scale ~ dpar(1.5,1)
prior.a <- prior.avg*prior.scale
prior.b <- (1-prior.avg)*prior.scale


# Addition to the likelihood model
p[i] ~ dbeta(prior.a[i],prior.b[i])
logit(prior.avg[i]) <- alpha + beta * Hunting[i]
prior.a[i] <- prior.avg[i]*prior.scale
prior.b[i] <- (1-prior.avg[i])*prior.scale

# Additional priors
prior.scale ~ dpar(1.5,1)
alpha ~ dunif(-20,20)
beta ~ dunif(-20,20)


#option (2)

mod2_string = " model {
     for (i in 1:length(y)) {
y[i] ~ dnorm(mu[z[i]], prec)
z[i] ~ dcat(omega)
}

#different means - two normal distributions
mu[1] ~ dnorm(-1.0, 1.0/100.0)
mu[2] ~ dnorm(1.0, 1.0/100.0) T(mu[1],) # ensures mu[1] < mu[2]

prec ~ dgamma(1.0/2.0, 1.0*1.0/2.0)
sig = sqrt(1.0/prec)

omega ~ ddirich(c(0.0, 1.0))
} "

