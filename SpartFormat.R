#melt data frame
library(reshape)
Obs <- melt(Obs,id=c("Site","Visit"))
names(Obs)[which(names(Obs)=="variable")] <- "Species"
names(Obs)[which(names(Obs)=="value")] <- "obs"

#Remove numbers from site and species names
Obs$Site <- as.numeric(gsub("Site","",Obs$Site))
Obs$Species <- as.numeric(gsub("Species","",Obs$Species))

#return data frame
return(Obs