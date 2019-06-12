getModelSummaries <- function(out){
  
  outCombined <- do.call(rbind,out)
  outSummaries <- ddply(outCombined,.(Type,Param),
                        summarise,
                        meanValue = mean(mean),
                        medValue = median(mean))
  
  return(outSummaries)
}
  