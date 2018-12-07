plotOccupancyHistory <- function(obs){
  require(ggplot)
  ggplot(obs)+
    geom_point(aes(x=Visit,y=Species,color=factor(obs)),size=2)+
    scale_colour_manual(values=c("black","white"))+
    theme(legend.position="none")+
    facet_wrap(~Site)
}
