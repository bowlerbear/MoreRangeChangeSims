getLQRepeatVisits <- function(Obs, NVisits, DetProb,
                            SiteDetEffects,YearDetEffects,IntDetEffects,
                            Scovariate,Tcovariate,Noise){

  
  #get rid of current Obs
  Occ <- Occ[,-grep("Visit",names(Occ))]
  
  #site and year random variation (observer variation)
  df <- unique(Obs[,c("Year","Site")])
  df$Noise <- rnorm(nrow(df),0,0.1)
  Obs$Noise <- df$Noise[match(interaction(Occ$Year,Occ$Site),
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

getLQRepeatAbundVisits <- function(Abund, NVisits, DetProb,
                                 SiteDetEffects,YearDetEffects,IntDetEffects,
                                 Scovariate,Tcovariate, Noise){
  
  #get rid of current obs
  Abund <- Abund[,-grep("Visit",names(Abund))]
  #names(Abund)[which(names(Abund)=="Occurence")] <- "Abund"
  
  #site and year random variation (observer variation)
  df <- unique(Abund[,c("Year","Site")])
  df$Noise <- rnorm(nrow(df),0,0.1)
  Abund$Noise <- df$Noise[match(interaction(Abund$Year,Abund$Site),
                                interaction(df$Year,df$Site))]
  
  #detection model (individual level) - linear predictor on logit scale
  lgtDetProb <- apply(Abund, 1, function(x) {
    logit(DetProb[x["Species"]]) +
      SiteDetEffects[x["Species"]] * ScovariateD[x["Site"]] + 
      YearDetEffects[x["Species"]] * Tcovariate[x["Year"]] +  
      IntDetEffects[x["Species"]] * ScovariateD[x["Site"]] * Tcovariate[x["Year"]]
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
  Abund$Occurence <- 0
  Abund$Occurence[Abund$Abund>0] <- 1
  
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


getLQSims <- function(x){
  
  if(samplingType=="Occurrence"){
    
    Obs <- getLQRepeatVisits(Occ=x, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate,Noise)
    return(Obs)
    
  } else if(samplingType=="Abundance"){
    
    Obs <- getLQRepeatAbundVisits(Abund=x, NVisits=NVisits, DetProb = DetProb, SiteDetEffects,YearDetEffects,IntDetEffects, Scovariate,Tcovariate,Noise)
    return(Obs)
  }
  
}
