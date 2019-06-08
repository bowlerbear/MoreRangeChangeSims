getJAGSformat <- function(Obs = Obs){
  
  #melt data frame
  require(reshape2)
  ObsMelted <- mel(Obs,id=c("Year","Species","Site","Occurence","DetProb"))
  names(ObsMelted)[which(names(ObsMelted)=="variable")]<-"Visit"
  names(ObsMelted)[which(names(ObsMelted)=="value")]<-"Obs"
  ObsMelted$Visit <- as.numeric(gsub("Visit","",ObsMelted$Visit))
  
  #calculate listlength per visit
  require(plyr)
  listlengthDF <- ddply(ObsMelted,.(Year,Site,Visit), summarise, nuSpecies=length(unique(Species[Obs>0])))
  
  #get species observation matrix
  obsMatrix <- acast(ObsMelted,Year+Species+Visit)
  
}