getRepeatVisits <- function(nVisits=5, PrObs = DetProb){
  
  Occ$NVisits <- rep(nVisits, nrow(Occ))
  
  #simulate the visits
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
  
  )
  
}