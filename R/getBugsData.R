getBugsData <- function(modelData){
  
  bugs.data <- list(nyear = length(unique(modelData$Year)),
                    nsite = length(unique(modelData$Site)),
                    year = modelData$Year,
                    site = modelData$Site,
                    y = modelData$Obs,
                    L = log(modelData$L)-median(log(modelData$L)),
                    nvisit = nrow(modelData),
                    sumX=sum(unique(modelData$Year)), 
                    sumX2=sum(unique(modelData$Year)^2))
  
  if("prevObs" %in% names(modelData)){
    bugs.data$prevObs <- modelData$prevObs
  }
  
  return(bugs.data)
  
}