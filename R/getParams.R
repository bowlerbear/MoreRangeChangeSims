getParams <- function(parameterDF,mysimNu){
  
  temp <- subset(parameterDF,SimNu==mysimNu) 
  
  #Number of sites in the landscape
  NSites <<- temp$Nsites
  
  #Number of species in the regional pool
  NSpecies <<- temp$Nspecies 
  
  #Number of time points
  NYears <<- temp$Nyears
  
  #Define average habitat quality of each site
  Scovariate <<- rnorm(NSites,0,temp$ScovariateVar)
  ScovariateD <<- Scovariate
  
  #Define time trend
  Tcovariate <<- (1:NYears)-1
  
  #sampling type
  samplingType <<- temp$samplingType
  
  #Define species responses
  speciesVariation <<- temp$speciesVariation
  
  if(speciesVariation==TRUE){
    OccProb <<- rnorm(NSpecies,temp$meanOccProb,temp$meanOccProb/5)#occupancy
    lambda <<- rnorm(NSpecies,temp$meanLambda,temp$meanLambda/10)#abundance
    Seffects <<- rnorm(NSpecies,temp$meanSeffects,temp$meanSeffects/5)##Define spatial quality effects
    Teffects <<- rnorm(NSpecies,temp$meanTeffects,temp$meanTeffects/5)#Effect of year
    Ieffects <<- rnorm(NSpecies,temp$meanIeffects,temp$meanIeffects/5)#Effect of year
  } else{
    OccProb <<- rep(temp$meanOccProb,NSpecies)
    lambda <<- rep(temp$meanLambda, NSpecies)
    Seffects <<- rep(temp$meanSeffects,NSpecies)
    Teffects <<- rep(temp$meanTeffects,NSpecies)
    Ieffects <<- rep(temp$meanIeffects,NSpecies)
  }
  
  #Define sampling effects
  
  #Number of visits to each site within each year
  NVisits <<- temp$Nvisits
  
  #are we subsetting to only positive results
  subsetPositive <<- temp$subsetPositive
  
  #Define average detection probability
  if(speciesVariation==TRUE){
    DetProb <<- rnorm(NSpecies,temp$meanDetProb,temp$meanDetProb/5) 
    SiteDetEffects <<- rnorm(NSpecies,temp$meanSiteDetEffects,temp$meanSiteDetEffects/5)
    YearDetEffects <<- rnorm(NSpecies,temp$meanYearDetEffects,temp$meanYearDetEffects/5)
    IntDetEffects <<- rnorm(NSpecies,temp$meanIntDetEffects,temp$meanIntDetEffects/5)
  }else{ 
    DetProb <<- rep(temp$meanDetProb,NSpecies)
    SiteDetEffects <<- rep(temp$meanSiteDetEffects,NSpecies)
    YearDetEffects <<- rep(temp$meanSiteDetEffects,NSpecies)
    IntDetEffects <<- rep(temp$meanSiteDetEffects,NSpecies)
  }
  
  #site/year detection noise (observer variation in effort)
  Noise <<- temp$Noise
  
  #Identify focal species
  focalSpecies <<- temp$focalSpecies
  
  #for focal species, set to fixed parameters
  OccProb[focalSpecies] <<- temp$meanOccProb
  lambda[focalSpecies] <<- temp$meanLambda
  Seffects[focalSpecies] <<- temp$meanSeffects
  Teffects[focalSpecies] <<- temp$meanTeffects
  Ieffects[focalSpecies] <<- temp$meanIeffects
  DetProb[focalSpecies] <<- temp$meanDetProb 
  SiteDetEffects[focalSpecies] <<- temp$meanSiteDetEffects
  YearDetEffects[focalSpecies] <<- temp$meanYearDetEffects
  IntDetEffects[focalSpecies] <<- temp$meanIntDetEffects
  
  
  #model type
  myModel <<- temp$myModel
  
}