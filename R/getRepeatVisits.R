# ####single time-slice#########################################
# #these functions take the true occurence history and assume sampling in a series of reapeated visits
# getRepeatVisits <- function(nVisits=NVisits, PrObs = DetProb){
#   
#   Occ$NVisits <- rep(nVisits, nrow(Occ))
#   
#   #simulate the visits
#   Obs <- sapply(rownames(Occ), function(site, Occ, PrObs = DetProb) {
#     nvisits <- Occ[site,"NVisits"]
#     occ <- unlist(Occ[site, names(Occ)!="NVisits"])
#     obs <- replicate(nvisits, rbinom(length(occ), 1, PrObs))
#     obs <- t(obs*occ)
#     colnames(obs) <- names(occ)
#     obs
#   }, Occ=Occ, simplify=FALSE)
#   
#   #unlist and reorganise the data frame    
#   Obs <- plyr::ldply(Obs)
#   names(Obs) <- gsub(".id", "Site", names(Obs))
#   
#   #Add Visit Number
#   Obs$Visit <- rep(1:nVisits,NSpecies)
#   
#   #melt data frame
#   library(reshape)
#   Obs <- melt(Obs,id=c("Site","Visit"))
#   names(Obs)[which(names(Obs)=="variable")] <- "Species"
#   names(Obs)[which(names(Obs)=="value")] <- "obs"
#   
#   #Remove numbers from site and species names
#   Obs$Site <- as.numeric(gsub("Site","",Obs$Site))
#   Obs$Species <- as.numeric(gsub("Species","",Obs$Species))
#   
#   
#   
#   #return data frame
#   return(Obs)
#   
# }
####multiple time slice##########################################################

getRepeatVisits <- function(Occ=Occ, NVisits=NVisits, DetProb=DetProb,
                            SiteDetEffects,YearDetEffects,IntDetEffects,
                            Scovariate,Tcovariate){
  
  #melt the array
  require(reshape2)
  Occ <- melt(Occ)
  names(Occ) <- c("Year","Species","Site","Occurence")
  
  #sort names
  Occ$Year <- as.numeric(gsub("Year","",Occ$Year))
  Occ$Species <- as.numeric(gsub("Species","",Occ$Species))
  Occ$Site <- as.numeric(gsub("Site","",Occ$Site))
  
  #detection model
  #linear predictor on logit scale
  lgtDetProb <- apply(Occ, 1, function(x) {
                              logit(DetProb) + 
                              SiteDetEffects[x["Species"]] * Scovariate[x["Site"]] + 
                              YearDetEffects[x["Species"]] * Tcovariate[x["Year"]] +  
                              IntDetEffects[x["Species"]] * Scovariate[x["Site"]] * Scovariate[x["Year"]]
                            })

  #convert into probabilities
  Occ$DetProb <- inv.logit(lgtDetProb)
  
  #null model
  #Occ$DetProb <- DetProb
  
  #apply detection probability to each visit - detection prob is constant within a year/site/species
  Obs <- replicate(NVisits,sapply(Occ$DetProb,function(x)rbinom(1,1,x)))
  
  #set to zero when the species was absent
  Obs[Occ$Occurence==0,] <- 0
  Obs <- as.data.frame(Obs)
  names(Obs) <- paste0("Visit",1:NVisits)
    
  #add repeat visits to the main data frame
  Occ <- as.data.frame(Occ)
  Occ <- cbind(Occ,Obs)
  
  #return data frame
  return(Occ)
 
}
