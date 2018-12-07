getRepeatVisits <- function(nVisits=5, PrObs = DetProb){
  
  if(length(PrOcc)==1) PrOcc <- rep(PrOcc, NPlots)
  if(length(PrObs)==1) PrObs <- rep(PrObs, NPlots)
  if(length(NVisits)==1) NVisits <- rep(NVisits, NPlots)
  
  Occ$NVisits <- rep(nVisits, nrow(Occ))
  
  #simulate the visits using binomial sampling
  Obs <- sapply(rownames(Occ), function(site, Occ, PrObs = DetProb) {
    nvisits <- Occ[site,"NVisits"]
    occ <- unlist(Occ[site, names(Occ)!="NVisits"])
    obs <- replicate(nvisits, rbinom(length(occ), 1, PrObs))
    obs <- t(obs*occ)
    colnames(obs) <- names(occ)
    obs
  }, Occ=Occ, simplify=FALSE)
  
  #unlist and reorganise the data frame    
  Obs <- plyr::ldply(Obs)
  names(Obs) <- gsub(".id", "Site", names(Obs))
  
  #Add Visit Number
  Obs$Visit <- rep(1:nVisits,NSpecies)

  #return data frame
  return(Obs)
  
}