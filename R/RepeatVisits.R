getRepeatVisits <- function(occ, nVisits=5, PrObs = DetProb){
  
  #extract metadata
  NPlots = nrow(occ)
  NSpecies = ncol(occ)
  
  #defense programming: number of visits per site
  if(length(nVisits)==1){
    Visits <- rep(nVisits, NPlots)
  }#allow each site to have a different number of visits
  else if(length(nVisits)==NPlots){
    Visits <- nVisits
  }else {
    stop("Number of visits and plots are not compatible")
  }
  
  #defense programming: detection probability per visit and site
  if(length(DetProb)==1) {
    PrOcc <- rep(DetProb, NPlots)
  }#allow each site to have a different number of visits
  else if(length(DetProb)==NPlots){
    PrOcc <- DetProb
  }else(length(NVisits)!=NPlots){
    stop("Number of visits and plot not compatible")
  }
  
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