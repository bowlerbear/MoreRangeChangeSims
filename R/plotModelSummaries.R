plotSummaries <- function(out){
  
  outCombined <- do.call(rbind,out)
  require(ggplot2)
  ggplot(subset(outCombined,Param %in% c("prop.p","prop.muZ","trend.muZ")))+
    geom_violin(aes(x=Type,y=mean),draw_quantiles=c(0.25,0.5,0.75))+
    facet_wrap(~Param)+
    theme_bw()+ylab("Estimate")
}
  