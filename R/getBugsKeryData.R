getBugsData <- function(modelData = modelData){
  
  bugs.data <- list(nyear = length(unique(modelData$Year)),
                    nsite = length(unique(modelData$Site)),
                    nabsvisit = length(unique(modelData$Visit)),
                    year = modelData$Year,
                    site = modelData$Site,
                    y = acast(modelData,Site~Year~Visit,value.var="Obs"),
                    L = log(modelData$L),
                    nvisit = nrow(modelData))
  
  return(bugs.data)
  
}