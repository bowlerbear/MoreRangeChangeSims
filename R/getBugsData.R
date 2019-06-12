getBugsData <- function(modelData = modelData){
  
  bugs.data <- list(nyear = length(unique(modelData$Year)),
                    nsite = length(unique(modelData$Site)),
                    Year = sort(unique(modelData$Year))-mean(unique(modelData$Year)),
                    year = modelData$Year,
                    site = modelData$Site,
                    y = modelData$Obs,
                    L = log(modelData$L)-mean(log(modelData$L)),
                    nvisit = nrow(modelData))
  
  return(bugs.data)
  
}