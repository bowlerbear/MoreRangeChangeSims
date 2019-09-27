# Format simulated recording data for use in sparta function
#
# This function takes the observation data generated under each scenario
# and carries out any formatting required to fit the dataset into the 
# sparta occupancy model function.

#melt data frame
getSpartaFormat<-function(Obs,focalSpecies,subsetPositive){
  
  #melt data frame
  require(reshape2)
  ObsMelted <- melt(Obs,id=names(Obs)[!grepl("Visit",names(Obs))])
  names(ObsMelted)[which(names(ObsMelted)=="variable")]<-"Visit"
  names(ObsMelted)[which(names(ObsMelted)=="value")]<-"Obs"
  ObsMelted$Visit <- as.numeric(gsub("Visit","",ObsMelted$Visit))
  ObsMelted$VisitID <- paste(ObsMelted$Year,ObsMelted$Site,ObsMelted$Visit,sep="_")
  
  #calculate listlength per visit
  require(plyr)
  listlengthDF <- ddply(ObsMelted,.(VisitID), summarise, nuSpecies=length(unique(Species[!is.na(Obs)&Obs>0])))
  ObsMelted$L <- listlengthDF$nuSpecies[match(ObsMelted$VisitID,listlengthDF$VisitID)]

  #subset to focal species
  ObsMeltedS <- subset(ObsMelted,Species==focalSpecies)
  
  #restrict to only occurences
  if(subsetPositive==T){
    ObsMeltedS <- subset(ObsMeltedS,L>0)
  }
  
  #return data frame 
  return(ObsMeltedS)

}