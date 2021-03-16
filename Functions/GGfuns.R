require(ggplot2)

# Correlations
GGcor <- function(x,y,name=NULL,annotations,labx,laby){
  require(ggalt)
  require(viridis)
  datatemp=data.frame(Obs=x,Pred=y,name=name)
  ggplot(datatemp, aes(Obs, Pred)) + 
    geom_point(shape=16, size=1, show.legend = FALSE) +
    stat_bkde2d(aes(fill=..level..), geom="polygon") +
    annotate(geom = "text", x =-Inf, y =Inf,
             label = paste("r=",round(cor(x,y,method="spearman"),2),sep=""),
             hjust =-0.1,vjust=1.3,size=10)+
    ggtitle(name)+
    scale_fill_viridis()+
    xlab(labx) +
    ylab(laby) +
    geom_abline(color="red")+
    theme(panel.background = element_rect(fill = 'grey85', colour = 'black'),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          axis.title=element_text(size=20),
          axis.text=element_blank(),axis.ticks=element_blank(),
          plot.title =element_text(size=30),
          legend.position = "none",
          plot.margin=unit(rep(-0.1,4), "cm"))
}

# Map predictions
GGmap=function(datRast,main,dir=1){
  dat=fortify(datRast)
  names(dat)=c("x","y","comp")
  ggplot()+coord_fixed(xlim=c(range(dat$x)[1],26), ylim=range(dat$y))+
    geom_tile(aes(x=x,y=y,fill=comp),data=dat)+
    scale_fill_distiller(name = "", palette ="Spectral" ,direction =dir)+
    ggtitle(main)+
    theme(panel.background = element_blank(),plot.title = element_text(size = 25, face = "bold"),
          legend.title = element_text(size = 15))
}