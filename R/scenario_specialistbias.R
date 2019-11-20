#get true occurences

getTrueSims <- function(){
  if(samplingType=="Occurrence"){
    Occ <- getOccuHistory(occProb = OccProb, Seffects, Teffects, Ieffects,
                          Tcovariate,Scovariate)
    return(Occ)
  } else if(samplingType=="Abundance"){
    Abund <- getAbundHistory(lambda = lambda, Seffects, Teffects, Ieffects,
                             Tcovariate,Scovariate)
    return(Abund)
  }
}

#get normal observer behaviour
getObs <- function(x){
  if(samplingType=="Occurrence"){
    Obs <- getRepeatVisits(Occ=x, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate,Noise)
    return(Obs)
    
  } else if(samplingType=="Abundance"){
    Obs <- getRepeatAbundVisits(Abund=x, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate,Noise)
    return(Obs)
  }
}

getbiasedSims <- function(propSites=0.5){
  
#run them and get simulation of a normal observer
temp <- getTrueSims()
normalCS <- getObs(temp)

#also get biased observer behaviour

#assume species detection probabbility vary with their occurence probability
speciesSpecialism <- rank(OccProb)-mean(rank(OccProb))
#species with wider distributions are assumed to be generalists

#the bias reduces their detection probability
DetProb <- inv.logit(logit(DetProb) - 0.05*speciesSpecialism)

#get observations of these people as well
biasedCS <- getObs(temp)

#a prop of people at the sites display this behaviour
sitesBiased <- sitesAutocorrelated <- sample(1:NSites,NSites*propSites)

#combine all data
Obs <- rbind(subset(normalCS,!Site %in% sitesBiased),
              subset(biasedCS, Site %in% sitesBiased))
             
return(Obs)

}

getbiaseddeclineSims <- function(biasEffect=0.15,propSites=1){
  
  #run them and get simulation of a normal observer
  temp <- getTrueSims()
  normalCS <- getObs(temp)
  
  #also get biased observer behaviour
  
  #assume some species are declining
  Teffects <- rnorm(NSpecies,0,0.04)
  
  #and this effects their detection probability
  speciesTrend<- rank(Teffects)-median(rank(Teffects))
  DetProb <- inv.logit(logit(DetProb) - biasEffect*speciesTrend)
  
  #get observations of these people as well
  biasedCS <- getObs(temp)
  
  #a prop of people at the sites display this behaviour
  sitesBiased <- sitesAutocorrelated <- sample(1:NSites,NSites*propSites)
  
  #combine all data
  Obs <- rbind(subset(normalCS,!Site %in% sitesBiased),
               subset(biasedCS, Site %in% sitesBiased))
  
  return(Obs)
  
}

