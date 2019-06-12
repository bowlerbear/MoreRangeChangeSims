#apply it
applyAutocorrelation <- function(Obs,prZero=1){
  Obs_auto <- Obs
  for(i in grep("Visit",names(Obs_auto))[-1]){
    Obs_auto[,i]<-ifelse((Obs_auto[,i]==1 & Obs_auto[,i-1]==1),rbinom(1,1,(1-prZero)),Obs_auto[,i])
  }
  return(Obs_auto)
}


#fix it by following 1s with NAs
fixAutocorrelation <- function(Obs_auto,addNAs="all"){
  #add NAs in the second visit after a detection
  Obs_fix <- Obs_auto
  for(i in grep("Visit",names(Obs_fix))[-1]){
    Obs_fix[,i]<-ifelse(Obs_fix[,i-1]==1|is.na(Obs_fix[,i-1]),NA,Obs_fix[,i])
  }
  return(Obs_fix)
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