getBugsData <- function(modelData = modelData){
  
  bugs.data <- list(nyear = length(unique(modelData$Year)),
                    nsite = length(unique(modelData$Site)),
                    year = modelData$Year,
                    site = modelData$Site,
                    y = modelData$Obs,
                    L = log(modelData$L),
                    nvisit = nrow(modelData))
  
  return(bugs.data)
  
}