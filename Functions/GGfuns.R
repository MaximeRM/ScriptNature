##%######################################################%##
#                                                          #
####                 PLOTTING FUNCTIONS                 ####
#                                                          #
##%######################################################%##


# Correlations
GGcor <- function(x, y, name=NULL, annotations, labx, laby) {
  datatemp <- data.frame(Obs=x, Pred=y, name=name)
  
  ggplot(datatemp, aes(Obs, Pred)) + 
    labs(title=name, x=labx, y=laby)
  geom_point(shape=16, size=1, show.legend = FALSE) +
    stat_bkde2d(aes(fill=..level..), geom="polygon") +
    annotate(geom = "text", x =-Inf, y =Inf,
      label = paste("r=", round(cor(x, y, method="spearman"), 2), sep=""),
      hjust = -0.1, vjust=1.3, size=10)+
    scale_fill_viridis()+
    geom_abline(color="red")+
    theme(panel.background = element_rect(fill = 'grey85', colour = 'black'),
      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
      axis.title = element_text(size=20),
      axis.text = element_blank(),axis.ticks=element_blank(),
      plot.title = element_text(size=30),
      legend.position = "none",
      plot.margin = unit(rep(-0.1,4), "cm"))
}

# Map predictions
GGmap <- function(datRast, main, dir=1){
  dat <- fortify(datRast)
  names(dat) <- c("x", "y", "comp")
  
  ggplot(dat)+
    labs(title=main)+
    coord_fixed(xlim=c(min(dat$x),26), ylim=range(dat$y))+
    geom_tile(aes(x=x, y=y, fill=comp))+
    scale_fill_distiller(name = "", palette ="Spectral" ,direction =dir)+
    theme(panel.background = element_blank(),plot.title = element_text(size = 25, face = "bold"),
      legend.title = element_text(size = 15))
}