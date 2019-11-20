#apply it
applyAutocorrelation <- function(Obs,propSites,autoProb){
  
  #list of sites visited by autocorrelated observers
  sitesAutocorrelated <- sample(1:NSites,NSites*propSites)
  
  Obs_auto <- Obs
  Obs_auto$autoSite <- ifelse(Obs_auto$Site %in% sitesAutocorrelated,1,0)
  #for each visit, find when it was seen twice in a row, 
  #and if it is a designed autocorrelation site
  #keep positive observations only with some probability
  for(i in grep("Visit",names(Obs_auto))[-1]){
    Obs_auto[,i]<-ifelse((Obs_auto[,i]==1 & 
                            Obs_auto[,i-1]==1 & 
                            Obs_auto$Site%in%sitesAutocorrelated),
                            rbinom(1,1,1-autoProb),
                            Obs_auto[,i])
  }
  return(Obs_auto)
}

applyFullAutocorrelation <- function(Obs,propSites){
  
  #list of sites visited by autocorrelated observers
  sitesAutocorrelated <- sample(1:NSites,NSites*propSites)
  
  Obs_auto <- Obs
  Obs_auto$autoSite <- ifelse(Obs_auto$Site %in% sitesAutocorrelated,1,0)
  
  #at these sites, set all Obs to 0 after the first 1
  for(i in 1:nrow(Obs_auto)){
    start = 0
    for(j in grep("Visit",names(Obs_auto))){
      if(start==0 & Obs_auto[i,j]==1 & Obs_auto$autoSite[i]==1){
        start = start + 1
      } else if(start>=1 & Obs_auto[i,j]==1 & Obs_auto$autoSite[i]==1){
        Obs_auto[i,j] = 0
      }
    }
  }
  
  Obs_auto <- Obs_auto[,-which(names(Obs_auto)=="autoSite")]
  
  return(Obs_auto)
}

fixAutocorrelation <- function(Obs_auto,bias="Autocorrelation"){
  
  modelData <- getSpartaFormat(Obs_auto,focalSpecies,subsetPositive)
  #add on detection on previous site
  modelData <- ddply(modelData,.(Year,Species,Site),function(x){
    x$prevObs <- 0
    for(i in 2:NVisits){
      x$prevObs[i-1] <- x$Ob[i-1]
      return(x)
    }
  })
  bugs.data <- getBugsData(modelData)
  myModel <- "basic_odm_fixedYear_LL_autocorrelation.txt"
  modelSummary <- runModel(bugs.data,modelData,focalSpecies,myModel)
  modelSummary$Bias <- bias
  
  return (modelSummary)
  
}


#fix it by following 1s with NAs
fixFullAutocorrelation <- function(Obs_auto,bias="FullAutocorrelation"){
  
  #add NAs in the second visit after a detection
  Obs_fix <- Obs_auto
  for(i in grep("Visit",names(Obs_fix))[-1]){
    Obs_fix[,i]<-ifelse(Obs_fix[,i-1]==1|is.na(Obs_fix[,i-1]),NA,Obs_fix[,i])
  }
  
  modelData <- getSpartaFormat(Obs_fix,focalSpecies,subsetPositive)
  #add on detection on previous site
  modelData <- ddply(modelData,.(Year,Species,Site),function(x){
    x$prevObs <- 0
    for(i in 2:NVisits){
      x$prevObs[i-1] <- x$Ob[i-1]
      return(x)
    }
  })
  bugs.data <- getBugsData(modelData)
  myModel <- "basic_odm_fixedYear_LL.txt"
  modelSummary <- runModel(bugs.data,modelData,focalSpecies,myModel)
  modelSummary$Bias <- bias
  
  return (modelSummary)
  
}

testAutocorrelation<- function(modelData){
  
  #create a null order
  modelData <- ddply(modelData,.(Year,Species,Site),function(x){
    x$nullOrder <- sample(x$Obs)
    return(x)
  })
  
  #compare number of consecutive observations
  consecutiveSummary <- ddply(modelData,.(Year,Species,Site),summarise,
                              nuNullConsecs = sum(diff(which(nullOrder==1))==1),
                              nuObsConsecs = sum(diff(which(Obs==1))==1))
  
  require(reshape2)
  ggplot(melt(consecutiveSummary[,c("nuNullConsecs","nuObsConsecs")]))+
    geom_violin(aes(x=variable,y=value))+scale_y_log10()+
    ylab("number of conseecutive observations")+xlab("Model")
  
  
}
