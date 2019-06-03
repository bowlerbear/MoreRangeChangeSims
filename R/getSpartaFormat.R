# Format simulated recording data for use in sparta function
#
# This function takes the observation data generated under each scenario
# and carries out any formatting required to fit the dataset into the 
# sparta occupancy model function.

#melt data frame
getSpartaFormat<-function(Obs){
library(reshape)
Obs <- melt(Obs,id=c("Site","Visit"))
names(Obs)[which(names(Obs)=="variable")] <- "Species"
names(Obs)[which(names(Obs)=="value")] <- "obs"

#return data frame
return(Obs)
}