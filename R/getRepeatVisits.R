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

getRepeatVisits <- function(Occ, NVisits, DetProb,
                            SiteDetEffects,YearDetEffects,IntDetEffects,
                            Scovariate,Tcovariate,Noise){
  
  #melt the array
  require(reshape2)
  Occ <- melt(Occ)
  names(Occ) <- c("Year","Species","Site","Occurence")
  
  #sort names
  Occ$Year <- as.numeric(gsub("Year","",Occ$Year))
  Occ$Species <- as.numeric(gsub("Species","",Occ$Species))
  Occ$Site <- as.numeric(gsub("Site","",Occ$Site))
  
  #site and year random variation (observer variation)
  df <- unique(Occ[,c("Year","Site")])
  df$Noise <- rnorm(nrow(df),0,0.01)
  Occ$Noise <- df$Noise[match(interaction(Occ$Year,Occ$Site),
                              interaction(df$Year,df$Site))]
  
  #detection model
  #linear predictor on logit scale
  lgtDetProb <- apply(Occ, 1, function(x) {
                              logit(DetProb[x["Species"]]) +
                              SiteDetEffects[x["Species"]] * Scovariate[x["Site"]] + 
                              YearDetEffects[x["Species"]] * Tcovariate[x["Year"]] +  
                              IntDetEffects[x["Species"]] * Scovariate[x["Site"]] * Tcovariate[x["Year"]]
                            })
  
  if(Noise==T){
    lgtDetProb <- lgtDetProb + Occ$Noise
  }
  
  #remove Noise variable from the dataframe
  Occ <- Occ[,-which(names(Occ)=="Noise")]
  
  
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


###########################################################################################

#assume abundance data
getRepeatAbundVisits <- function(Abund, NVisits, DetProb,
                            SiteDetEffects,YearDetEffects,IntDetEffects,
                            Scovariate,Tcovariate, Noise){
  
  #melt the array
  require(reshape2)
  Abund <- melt(Abund)
  names(Abund) <- c("Year","Species","Site","Abund")
  
  #sort names
  Abund$Year <- as.numeric(gsub("Year","",Abund$Year))
  Abund$Species <- as.numeric(gsub("Species","",Abund$Species))
  Abund$Site <- as.numeric(gsub("Site","",Abund$Site))
  
  #site and year random variation (observer variation)
  df <- unique(Abund[,c("Year","Site")])
  df$Noise <- rnorm(nrow(df),0,0.01)
  Abund$Noise <- df$Noise[match(interaction(Abund$Year,Abund$Site),
                              interaction(df$Year,df$Site))]
  
  #detection model (individual level) - linear predictor on logit scale
  lgtDetProb <- apply(Abund, 1, function(x) {
      logit(DetProb[x["Species"]]) +
      SiteDetEffects[x["Species"]] * Scovariate[x["Site"]] + 
      YearDetEffects[x["Species"]] * Tcovariate[x["Year"]] +  
      IntDetEffects[x["Species"]] * Scovariate[x["Site"]] * Tcovariate[x["Year"]]
  })
  
  if(Noise==T){
    lgtDetProb <- lgtDetProb + Abund$Noise
  }
  
  #remove Noise variable from the dataframe
  Abund <- Abund[,-which(names(Abund)=="Noise")]
  
  #convert into probabilities
  Abund$DetProb <- inv.logit(lgtDetProb)
  
  #null model
  #Occ$DetProb <- DetProb
  
  #convert to detection probability of species (see at least individual)
  Abund$DetProb <- apply(Abund,1,function(x)(1-((1-x["DetProb"])^x["Abund"])))
  
  #apply detection probability to each visit - detection prob is constant within a year/site/species
  Obs <- replicate(NVisits,sapply(Abund$DetProb,function(x)rbinom(1,1,x)))
  
  #convert abundance to occurence
  names(Abund)[which(names(Abund)=="Abund")] <- "Occurence"
  Abund$Occurence[Abund$Occurence>0] <- 1
  
  #add on replciate visits
  Obs[Abund$Occurence==0,] <- 0#not necessary
  Obs <- as.data.frame(Obs)
  names(Obs) <- paste0("Visit",1:NVisits)
  
  #add repeat visits to the main data frame
  Abund <- as.data.frame(Abund)
  Abund <- cbind(Abund,Obs)
  
  #return data frame
  return(Abund)
  
}




###########################################################################################
