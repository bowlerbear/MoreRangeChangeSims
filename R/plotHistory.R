plotHistory <- function(obs){
  require(ggplot2)
  ggplot(obs)+
    geom_point(aes(x=Visit,y=Species,color=factor(obs)),size=2)+
    scale_colour_manual(values=c("white","black"))+
    theme(legend.position="none")+
    facet_wrap(~Site)
}
