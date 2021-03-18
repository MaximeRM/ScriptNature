# Required packages
libs <- c("ade4","ape","SCGLR","ggplot2","raster","maptree","ggalt","viridis")
for(l in libs) {
  library(l, character.only = TRUE)
}

# Required functions
source("Functions/TAG.R")
source("Functions/GGfuns.R")

#################################################
# LOAD & ORGANIZE DATA 

# Abundance data (available upon request see XXXXXXXXXXXX)
abond <- readRDS("Data/CoForData.rds")

# Species names
nY <- names(abond)[2:194]
# Climate predictors
nX <- names(abond)[203:226]
# Additional variables
nT <- c("HWSD","Anthr2", "I(log(Anthr2))")


#################################################
# SCGLR analyses

# Formula
form <- multivariateFormula(nY, c(nX), A = nT)
# Type of family distribution
family <- rep("poisson", length(nY))

# Crossvalidation to select the number of climate components: Time consuming!
cv <- scglrCrossVal(form, data=abond,K=6, offset=abond$offset, nfolds=10,
                    family=family, crit=list(tol=1e-6,maxit=1000),
                    method=methodSR(l=1,s=0.1,bailout = 1000))
mean.crit <- colMeans(log(cv))
kOpt <- which.min(mean.crit)-1 # kOpt=3 in the present study

# Run scglr on the whole dataset
genus.scglr <-scglr(form, data = abond,K = kOpt, offset = abond$offset,
                    family = family, crit = list(tol = 1e-6, maxit = 1000),
                    method = methodSR(l = 1,s = 0.1,bailout = 1000,maxiter = 100))
plot(genus.scglr,threshold=0.1,labels.size=0.6,labels.auto=FALSE,plane=c(1,2))

### Assess prediction error with spatial cross validations: Time consuming!

# Create spatial clusters
ngroup <- 40 # Number of spatial clusters, i.e. number of cal/val datasets
tmp_clust <- abond[,c("lon", "lat", "id_cell")]
mdist <- dist(tmp_clust[c("lon","lat")])
hc <- hclust(mdist, method="ward.D")
tmp_clust$Clust_100km = cutree(hc, k=ngroup)

# plot clusters
plot(tmp_clust$lon, tmp_clust$lat, col=sample(rainbow(ngroup))[as.factor(tmp_clust$Clust_100km)],
      pch=20,xlab="lon",ylab="lat",asp=1,cex=0.5)

# Predictions on the validation dataset
repet <- ngroup
gpeval <- tmp_clust$Clust_100km
plots <- 1:nrow(abond)

mat_pred_val = as.matrix(abond[, nY])
for (i in 1:repet) {
  message(i, "/", repet)
  plots_val <- plots[gpeval == i]
  plots_cal <- plots[gpeval != i]
  genus.scglr2 <- scglr(form,data=abond,K=kOpt,offset=abond$offset,subset=plots_cal,
                        family=family,crit=list(tol=1e-6,maxit=1000),
                        method=methodSR(l=4,s=0.25,bailout = 1000,maxiter = 100,epsilon = 1e-9))
  ## Validation matrix
  x_new <-  model.matrix(form, data = abond[plots_val,], rhs = 1:length(form)[2])
  ## Predicting abundances on validation dataset
  mat_pred_val[plots_val,] <- multivariatePredictGlm(x_new,family=rep("poisson",length(nY)),
                                                     beta=as.matrix(genus.scglr2$beta),
                                                     offset=abond$offset[plots_val])
}

# Summary statistics at taxa level
spear <- diag(cor(as.matrix(abond[, nY]), mat_pred_val,method="spearman"))
summary(spear)

# Summary statistics at community level
dudiCA <- dudi.coa(abond[,nY], scann =FALSE, nf = 10) # Correspondence analysis (CA) done on observations
compowholeCA <- suprow(dudiCA, mat_pred_val[, nY]) # Projecting predictions in the CA performed on observations 
diag(cor(dudiCA$li, compowholeCA$lisup,method="spearman"))


#################################################
# Predicting abundances over the study area using the fitted SCGLR model

# Load study area with predictors
GridTot <- readRDS("Data/GridTot.rds")


# Predictions over the full study area
xTot <- model.matrix(form, data = GridTot, rhs = 1:length(form)[2])
prediction <- data.frame(multivariatePredictGlm( xTot,family = rep("poisson", length(nY)),beta = as.matrix(genus.scglr$beta)))
prediction <- data.frame(prediction, lon = GridTot$lon,lat = GridTot$lat)

# Set to zero all prediction close to zero (taking the minimum observed abundance as threshold) to avoid deleterious effect on CA
abondPerHa <- abond[,nY]/abond$offset
Ab <- apply(abondPerHa,2,function(x) min(x[x>0]))
thresholdAb <- min(Ab)
prediction[,nY][prediction[,nY] < thresholdAb] <- 0

# Floristic analysis 
datFlo <- prediction[,nY]
dudiPred <- dudi.coa(datFlo, scannf=FALSE, nf=5)
dataCA <- data.frame(prediction[, c("lon","lat")],dudiPred$li[,1:3])
dataCA <- dataCA[GridTot$CCI != 160,]
stackCA <- rasterFromXYZ(dataCA)
plot(stackCA)

writeRaster(stackCA,file="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/Round3/Soumission/Revision/Submission/Dataverse/Floristic_CA_scores/FloristicCAscores.tif",
            format="GTiff",overwrite=TRUE)

# Map functional traits
source("Functions/TAG.R")
tabTrait <- readRDS("Data/traitTab.rds")
cwmTrait  <- as.data.frame(TAG(tabTrait[nY,],prediction[,nY])$tag)
dataTrait <- data.frame(prediction[, c("lon","lat")],cwmTrait)
dataTrait <- dataTrait[GridTot$CCI != 160,]
stackTrait <- rasterFromXYZ(dataTrait)
plot(stackTrait)

writeRaster(stackTrait,file="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/Round3/Soumission/Revision/Submission/Dataverse/Functional_composition/FunctionalComposition.tif",
            format="GTiff",overwrite=TRUE)


# Build floristic types
Naxes <- 5
ngroup <- 10
TreeClust <- hclust(dist(dudiPred$li[,1:Naxes]),method = "ward.D")
SubTree <- clip.clust(TreeClust,k = ngroup,data=prediction)
extract <- function(cl, m){
  res <- data.frame(prediction[as.integer(m), c("lon","lat")], clFac=cl)
}
classPredict <- mapply(extract, names(SubTree$membership), SubTree$membership, SIMPLIFY = FALSE)
classPredict <- do.call(rbind, classPredict)
filt <- GridTot$CCI != 160
classPredict <- classPredict[paste(classPredict$lon,classPredict$lat)%in%paste(GridTot$lon[filt],GridTot$lat[filt]),]
typeFlo <- rasterFromXYZ(classPredict)
plot(typeFlo,col=rainbow(10)[sample(1:10)])

writeRaster(typeFlo,file="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/Round3/Soumission/Revision/Submission/Dataverse/Floristic_types/FloristicType.tif",
            format="GTiff",overwrite=TRUE)

pred <- prediction[GridTot$CCI != 160,]
write.csv(pred,file="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/Round3/Soumission/Revision/Submission/Dataverse/Predicted_abundances/PredictedAbundance.csv",
          row.names = F)


#################################################
# Predicting vulnerability

##########################
# Load present and future climate components
CompAfricaPresRast <- stack("Data/CompAfricaPresent.tif")



# #############################################################################
## SENSITIVITY CLIMATE

#Function sd weighted (adapted from SDMTools package)
wt.sd <- function(x,wt) {
    s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
    wtmean= sum(wt * x)/sum(wt) #return the mean
    wtvar=sum(wt *(x-wtmean)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2)))
    return( sqrt(wtvar) ) #return the standard deviation
} 

# Add components to observations
CompObs <- as.data.frame(raster::extract(CompAfricaPresRast,abond[,c("lon","lat")]))
names(CompObs) <- c("CC1","CC2", "CC3")
ObsTaxa=cbind(abond[,nY],CompObs)
  
CC1obssd=c();CC2obssd=c();CC3obssd=c()
for (i in 1:length(nY)){
    # Observations
    tempObs=ObsTaxa[,nY[i]]
    CC1obssd[i]=wt.sd(ObsTaxa$CC1,tempObs)
    CC2obssd[i]=wt.sd(ObsTaxa$CC2,tempObs)
    CC3obssd[i]=wt.sd(ObsTaxa$CC3,tempObs)
}
datSensitivitySPobs <- data.frame(CC1obssd,CC2obssd,CC3obssd,row.names = nY)
datSensitivitySPobs$meanCCsd <- apply(datSensitivitySPobs,1,mean)

# Mean sensitivity at the community level
TAGres <- as.data.frame(TAG(datSensitivitySPobs,PredTaxa[,nY])$tag)
MeanSensitivityClim  <- rasterFromXYZ(cbind(PredTaxa[,c("lon","lat")],1/TAGres$meanCCsd))
plot(MeanSensitivityClim)
  

# #############################################################################
## EXPOSURE TO CLIMATE CHANGE

# Function to compute dist between datasets
FUNDIST=function(x){
  if(any(is.na(x))){
    res=NA
  }else{
    res=dist(rbind(x[1:3],x[4:6]))
  }
  return(res)
}

CompAfricaPresRastCrop=crop(CompAfricaPresRast,extent(MeanSensitivityClim))
coordCrop=coordinates(CompAfricaPresRastCrop)
  
FunDistPixel=function(fileFut){
  dataClimFut=stack(fileFut)
  dataClimFut=crop(dataClimFut,extent(MeanSensitivityClim))
  names(dataClimFut)=c(paste0("C",1:19),c("meanET0","meanCWB","sumCWD","maxCWD","MCWD"))
  dataClimFutDat=as.data.frame(dataClimFut)
  # Futur component
  scaledVar=as.matrix(sweep(dataClimFutDat,2,meanPres,FUN="-"))%*%sdPres
  CompFutur=as.data.frame(as.matrix(scaledVar)%*%as.matrix(ResList$genus.scglr$u))
  CompFuturRast=rasterFromXYZ(cbind(coordCrop,CompFutur))
  # Euclidean distances between present and future climatic conditions
  Pres=as.data.frame(CompAfricaPresRastCrop)
  Fut=as.data.frame(CompFuturRast)
  FutMinusPres=Fut-Pres
  aa=cbind(Pres,Fut)
  distClim=apply(aa,1,FUNDIST)
  CompAfricaDistRast=CompAfricaPresRastCrop[[1]]
  CompAfricaDistRast[]=distClim
  return(list(CompAfricaDistRast,FutMinusPres,dataClimFutDat))
}

###### LOOP ON INDIVIDUAL SCENARIOS
listClimateModel=read.csv(paste0(PathClim,"listClimateModel.csv"))
  
resRCP45_2085=stack()
for (i in 1:nrow(listClimateModel)){
  fileFut45_85=paste0(PathClim,"dataClimFut_",listClimateModel[i,2],"_",listClimateModel[i,1],"_RCP45_2085.grd")
  #
  temp=FunDistPixel(fileFut45_85)
  resRCP45_2085=stack(resRCP45_2085,temp[[1]])
}
names(resRCP45_2085)=paste(toupper(sapply(strsplit(as.character(listClimateModel$cc),split="_"),"[",1)),c("",1:4,1,1:2,1:10))

#Calculate mean giving equal weights to the 18 models
r_mean45_2085 <- mean(resRCP45_2085)
MeanExpositionClim <- resample(r_mean45_2085,MeanSensitivityClim)
MeanExpositionClim <- crop(MeanExpositionClim,MeanSensitivityClim)
MeanExpositionClim[is.na(MeanSensitivityClim)]=NA
plot(MeanExpositionClim)
  
# #############################################################################
## ADAPTATION
library(picante)
library(entropart)
library(adephylo)
library(phylobase)
  
# List of Genera in dataset
Gendata <- sapply(strsplit(DataList$nY,split="_"),"[",1)
  
# Phylogenetic tree: Janssens, S. B. et al (2020). A large-scale species level dated angiosperm phylogeny for evolutionary and ecological analyses. Biodiversity data journal, 8.
PhyloTree <- read.nexus("Data/bdj-08-e39677-s005.tre")
taxa <- PhyloTree$tip.label
taxaGen <- sapply(strsplit(taxa,split="_"),"[",1)

# Prune the tree with only the genera of the dataset
filtsp=taxa%in%nY
TipToKeepSp=taxa[filtsp][!duplicated(taxaGen[filtsp])]
filtgen=taxaGen%in%Gendata & !filtsp & !taxaGen%in%taxaGen[taxa%in%nY]
TipToKeepGen=taxa[filtgen][!duplicated(taxaGen[filtgen])]
Tree=keep.tip(PhyloTree,tip=c(TipToKeepSp,TipToKeepGen))

  # Change names to genus level
  Tree$tip.label=sapply(strsplit(Tree$tip.label,split="_"),"[",1)
  length(Tree$tip.label)
  
  # Change names in Prediction
  dataPred=PredTaxa[,nY]
  dataObsOri=as.data.frame(DataList$abondsub)
  dataObs=dataObsOri[,nY]
  # Sum individuals per genus
  dataPredt=as.data.frame(t(dataPred))
  dataPredt$genus=Gendata
  dataPredt=data.table(dataPredt)
  dataPredtg=dataPredt[,lapply(.SD, sum, na.rm=TRUE),by=genus]
  genus=dataPredtg$genus
  dataPredtg$genus=NULL
  dataPredg=as.data.frame(t(dataPredtg))
  colnames(dataPredg)=genus
  
  dataObst=as.data.frame(t(dataObs))
  dataObst$genus=Gendata
  dataObst=data.table(dataObst)
  dataObstg=dataObst[,lapply(.SD, sum, na.rm=TRUE),by=genus]
  genus=dataObstg$genus
  dataObstg$genus=NULL
  dataObsg=as.data.frame(t(dataObstg))
  colnames(dataObsg)=genus
  
  dataPred=dataPredg[,Tree$tip.label] 
  dataObs=dataObsg[,Tree$tip.label] 
  dim(dataPred)
  dim(dataObs)
  
  # Test phylogenetic signal in component effects
  # Effect of component on each taxa
  taxScore=ResList$genus.scglr$gamma
  taxScoret=as.data.frame(t(taxScore))
  taxScoret$genus=Gendata
  taxScoret=data.table(taxScoret)
  taxScoretg=taxScoret[,lapply(.SD, mean, na.rm=TRUE),by=genus]
  genus=taxScoretg$genus
  taxScoretg$genus=NULL
  taxScoretg=as.data.frame(taxScoretg)
  rownames(taxScoretg)=genus
  
  taxScore=taxScoretg[Tree$tip.label,2:4]
  dim(taxScore)
  
  tree_in_phylo4 <- as(Tree,"phylo4")
  aa=proxTips(tree_in_phylo4,method="oriAbouheif")
  abouheif.moran(taxScore, aa,nrepet = 9999)
  
  # Chao index with q=1
  # dataPredChao=t(apply(dataPred,1,function(x) as.ProbaVector(x)))
  # ChaoDiv=apply(dataPredChao,1,function(x) ChaoPD(x,q=1,Tree))
  # save(ChaoDiv,file="ChaoDiv")
  # dataObsChao=t(apply(dataObs,1,function(x) as.ProbaVector(x)))
  # ChaoDivObs=apply(dataObsChao,1,function(x) ChaoPD(x,q=1,Tree))
  # save(ChaoDivObs,file="ChaoDivObs")
  load("ChaoDiv")
  load("ChaoDivObs")
  ChaoDivRast=rasterFromXYZ(cbind(PredTaxa[,c("lon","lat")],ChaoDiv))
  ChaoDivRastObs=rasterFromXYZ(cbind(dataObsOri[,c("lon","lat")],ChaoDivObs))
  ChaoPlot=GGfun(ChaoDivRast,main="Adaptation",sizeMain = 5,dir=1,stretch = "lin")
  ggsave("AdaptationChao.png",ChaoPlot,width=14,height=8)
  
  Adaptation=ChaoDivRast
  writeRaster(Adaptation,"Guillaume/Adaptation.grd",overwrite=T)
  
  # Validation
  aa=crop(ChaoDivRast,ChaoDivRastObs)
  bb=stack(aa,ChaoDivRastObs)
  cc=na.omit(as.data.frame(bb))
  ChaoplotCor=GGcor(cc$ChaoDivObs,cc$ChaoDiv,name="Chao",labx="Observed",laby="Predicted")
}

# #############################################################################
## VULNERABILITY 

########################################################
# Normalize the three components of vulnerability
MeanSensitivityNorm=calc(MeanSensitivityClim,scales::rescale)
MeanExpositionTotNorm45_2055=calc(MeanExpositionClimStack$r_mean45_2055,scales::rescale)
MeanExpositionTotNorm85_2055=calc(MeanExpositionClimStack$r_mean85_2055,scales::rescale)
MeanExpositionTotNorm45_2085=calc(MeanExpositionClimStack$r_mean45_2085,scales::rescale)
MeanExpositionTotNorm85_2085=calc(MeanExpositionClimStack$r_mean85_2085,scales::rescale)
AdaptationNorm=calc(Adaptation,scales::rescale)

# Calculate vulnerability
VulnerabilityClim45_2055=MeanSensitivityNorm+MeanExpositionTotNorm45_2055-AdaptationNorm
VulnerabilityClim85_2055=MeanSensitivityNorm+MeanExpositionTotNorm85_2055-AdaptationNorm
VulnerabilityClim45_2085=MeanSensitivityNorm+MeanExpositionTotNorm45_2085-AdaptationNorm
VulnerabilityClim85_2085=MeanSensitivityNorm+MeanExpositionTotNorm85_2085-AdaptationNorm

VulnerabilityMaps=stack(MeanSensitivityClim,
                        MeanExpositionClimStack,
                        Adaptation,
                        VulnerabilityClim45_2055,
                        VulnerabilityClim85_2055,
                        VulnerabilityClim45_2085,
                        VulnerabilityClim85_2085)

writeRaster(VulnerabilityMaps,filename="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/AccesData/Dataverse/Climate_vulnerability/ClimateVulnerability.tif",format="GTiff",overwrite=TRUE)

#####################################################
# Figure vulnerability

MainSubPlotSize=9

# GCO pourquoi pas GGfunPaper comme les autres ?
# stretch c'est plutot cut non ?
VulnerabilityClimNormplot45_2055=GGfun(VulnerabilityClim45_2055,sizeMain =16,dir=-1,stretch = "lin")+
  labs(title = "A Vulnerability to climate change")
VulnerabilityClimNormplot85_2055=GGfun(VulnerabilityClim85_2055,sizeMain =16,dir=-1,stretch = "lin")+
  labs(title = "A Vulnerability to climate change")
VulnerabilityClimNormplot45_2085=GGfun(VulnerabilityClim45_2085,sizeMain =16,dir=-1,stretch = "lin")+
  labs(title = "A Vulnerability to climate change")
VulnerabilityClimNormplot85_2085=GGfun(VulnerabilityClim85_2085,sizeMain =16,dir=-1,stretch = "lin")+
  labs(title = "A Vulnerability to climate change")

MeanSensitivityTOTplotNorm=GGfunPaper(MeanSensitivityNorm,
                                      Title ="C Sensitivity" ,
                                      SubTitle ="(to current climate gradients)",
                                      sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")

MeanMeanExpositionTotPlotNorm45_2055=GGfunPaper(MeanExpositionTotNorm45_2055,
                                                Title ="D Exposure",
                                                SubTitle ="(to forecasted climate change)",
                                                sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")
MeanMeanExpositionTotPlotNorm85_2055=GGfunPaper(MeanExpositionTotNorm85_2055,
                                                Title ="D Exposure",
                                                SubTitle ="(to forecasted climate change)",
                                                sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")
MeanMeanExpositionTotPlotNorm45_2085=GGfunPaper(MeanExpositionTotNorm45_2085,
                                                Title ="D Exposure",
                                                SubTitle ="(to forecasted climate change)",
                                                sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")
MeanMeanExpositionTotPlotNorm85_2085=GGfunPaper(MeanExpositionTotNorm85_2085,
                                                Title ="D Exposure",
                                                SubTitle ="(to forecasted climate change)",
                                                sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")


MeanAdaptationTOTplotNorm=GGfunPaper(AdaptationNorm,
                                     Title ="E Adaptive capacity",
                                     SubTitle ="(phylogenetic diversity)",
                                     sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")

# FigFinal={
#   (VulnerabilityClimNormplot)/(MeanSensitivityTOTplotNorm+MeanMeanExpositionTotPlotNorm+MeanAdaptationTOTplotNorm)+
#     plot_layout(nrow =2,height = c(2, 1.5))
# }
# ggsave("FigVulnerabilityClim.png",FigFinal,width=10,height=6.5)

############################################"
#### Anthropic vulnerability
############################################"
#### Anthropic vulnerability
Anthrpresent=raster("/home/rejou/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/OLD/FinalAnthropoIndex/AnthropoIndex2000__3Jan.tif")

Anthrpresentnew=raster("/home/rejou/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/OLD/FinalAnthropoIndex/AnthropoIndex2000_28_Apr2020.tif")
Anthr2085new=raster("/home/rejou/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/OLD/FinalAnthropoIndex/AnthropoIndex2085_28_Apr2020.tif")

ext=extent(Anthrpresentnew)
ext@xmin=8
AnthrMaps=stack(Anthrpresentnew,Anthr2085new)
names(AnthrMaps)=c("Anthr_present",
                   "Anthr_2085")
AnthrMaps=crop(AnthrMaps,ext)
plot(AnthrMaps)
writeRaster(AnthrMaps,"~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/AccesData/Dataverse/Anthropogenic_pressure/Anthropogenicpressure.tif",
            format="GTiff",overwrite=TRUE)


Anthrpresentnew=crop(Anthrpresentnew,Anthrpresent)
Anthr2085new=crop(Anthr2085new,Anthrpresent)

presAnthr=resample(Anthrpresentnew,MeanSensitivityClim)
presAnthr=crop(presAnthr,MeanSensitivityClim)
presAnthr[is.na(MeanSensitivityClim)]=NA

VulnerabilityAnthr=resample(Anthr2085new,MeanSensitivityClim)
VulnerabilityAnthr=crop(VulnerabilityAnthr,MeanSensitivityClim)
VulnerabilityAnthr[is.na(MeanSensitivityClim)]=NA

rangeAnthr=range(c(range(presAnthr[],na.rm = T),range(VulnerabilityAnthr[],na.rm = T)))

present=GGfun(presAnthr,sizeMain =16,dir=-1,stretch = "none",samelegend = T,rangeleg =rangeAnthr,lowhighVal=T)+
  labs(title = "A Current anthropogenic pressure")
futur=GGfun(VulnerabilityAnthr,sizeMain =16,dir=-1,stretch = "lin",samelegend = T,rangeleg =rangeAnthr,lowhighVal=T)+
  labs(title = "B Anthropogenic pressure projected in 2085")
pp=present+futur
ggsave("CurrentAndFutureAnthr.pdf",pp,width=10,height=5)

VulnerabilityAnthrplot=GGfunPaper(VulnerabilityAnthr,
                                  Title="",
                                  SubTitle = "",
                                  sizeMain =10,dir=-1,stretch = "lin")
VulnerabilityAnthrplot=GGfun(VulnerabilityAnthr,sizeMain =MainSubPlotSize,dir=-1,stretch = "lin")+
  labs(title =NULL)


###################################################
# LINK BETWEEN VULNERABILITY AND FUNCTIONS
CWD=stack("~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/AccesData/Dataverse/Functional_composition/FunctionalComposition.tif")
CWD=resample(CWD,MeanSensitivityClim)
CWD=crop(CWD,MeanSensitivityClim)
CWD[is.na(MeanSensitivityClim)]=NA

aa=na.omit(as.data.frame(stack(CWD,VulnerabilityAnthr,VulnerabilityClim45_2085)))
names(aa)=c("WD","Deciduousness","diamMax","Anthr","Clim")
str(aa)
bb=as.data.frame(cor(aa,method="spearman"))
write.csv(bb,"CorrelationVulnerabilityTrait.csv")

##################################################
# GLOBAL VULNERABILITY FIGURE 3
source("Guillaume/DataAndFuns/Myfunctions.R")
source("Guillaume/DataAndFuns/patch.R")
library(colorplaner)# remotes::install_github("wmurphyrd/colorplaner") # the 23/04/2020

pp=stack(VulnerabilityClim45_2085,VulnerabilityAnthr)
bb=as.data.frame(pp)
cor.test(bb$layer,bb$AnthropoIndex2085_28_Apr2020)




VulnerabilityAnthrNorm=calc(VulnerabilityAnthr,scales::rescale)
VulnerabilityClimNorm45_2055=calc(VulnerabilityClim45_2055,scales::rescale)
VulnerabilityClimNorm85_2055=calc(VulnerabilityClim85_2055,scales::rescale)
VulnerabilityClimNorm45_2085=calc(VulnerabilityClim45_2085,scales::rescale)
VulnerabilityClimNorm85_2085=calc(VulnerabilityClim85_2085,scales::rescale)

aa45_2055=stack(VulnerabilityClimNorm45_2055,VulnerabilityAnthrNorm)
aa85_2055=stack(VulnerabilityClimNorm85_2055,VulnerabilityAnthrNorm)
aa45_2085=stack(VulnerabilityClimNorm45_2085,VulnerabilityAnthrNorm)
aa85_2085=stack(VulnerabilityClimNorm85_2085,VulnerabilityAnthrNorm)

globalplot=function(aa,MeanMeanExp,MainTitle){
  dat=as.data.frame(aa,xy=T)
  names(dat)=c("x","y","VulnerabilityClim","VulnerabilityAnthr")
  # Apply a strech to variables
  dat$VulnerabilityClim=Strech(dat$VulnerabilityClim)
  dat$VulnerabilityAnthr=Strech(dat$VulnerabilityAnthr)
  
  paperTheme <- theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                      panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      axis.title=element_blank(),
                      plot.subtitle = element_text(hjust = 0.5,size=6),
                      axis.line=element_blank(),
                      axis.text.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      legend.position="none",
                      panel.border=element_blank())
  
  # Start Figure
  VulnTot <- ggplot()+
    labs(
      title=MainTitle,
      x="long", y="lat"
    )+
    
    # mise en forme
    paperTheme+
    theme(panel.background = element_blank(),plot.title = element_text(size = 12, face = "bold",hjust = 0.5),
          legend.title = element_text(size = 12),legend.position = "right")+
    
    # coordonnées
    # coord_fixed(xlim=c(min(dat$x),26), ylim=range(dat$y))+
    coord_fixed(xlim=c(8.1,27), ylim=c(-6,5))+
    
    # vulnerabilité
    geom_tile(aes(x = x, y = y, fill = VulnerabilityClim, fill2 = VulnerabilityAnthr),data=dat)+
    
    # domaine de calibration
    geom_tile(aes(x=lon,y=lat),fill="grey60",data=outofcal)+
    
    # échelle de couleur
    scale_fill_colorplane(name="toto",axis_title ="Vulnerability to \nclimate change", 
                          axis_title_y = "B Anthropogenic \npressure in 2085", 
                          color_projection = YUV_projection, Y=0.65,
                          na.color = NA,
                          breaks =c(0.1,0.9), breaks_y =c(0.1,0.9), 
                          labels =c("Low","High"), labels_y =c("Low","High"))+
    
    # légende
    guides(fill=guide_colorplane(
      title.theme =element_text(size=30,color="white"),
      label.theme = element_text(colour = "black",size=7),
      label.y.theme = element_text(colour = "black",size=7,angle=90),
      axis.title.y.position = "right",
      title.vjust =2.5,
      axis.title.theme = element_text(size=9,face="bold"),
      axis.title.y.theme = element_text(size=9,angle=90,face="bold")
    ))+
    
    # contour des pays
    # geom_sf(data=afrique,color="black",fill=NA)
    geom_path(aes(x=long,y=lat,group=group),color="black",data=country.point)
  
  # pour débugger ajoute un fond rose pour voir ce qui se passe
  if(FALSE) {
    debug_thm <- theme(plot.background = element_rect(fill="pink",color="black"))
  } else {
    debug_thm <- theme()
  }
  
  Ylim=ylim(-6,5)
  
  A <- VulnTot+theme(legend.position = "none")+Ylim
  LA <- cowplot:::get_legend(VulnTot)
  
  B <- VulnerabilityAnthrplot+paperTheme+Ylim+theme(legend.position = "none")
  LB <- cowplot:::get_legend(VulnerabilityAnthrplot+theme(legend.position = "right")) 
  
  C <- MeanSensitivityTOTplotNorm+Ylim
  D <- MeanMeanExp+Ylim
  E <- MeanAdaptationTOTplotNorm+Ylim
  
  plus <- grid::textGrob("+",gp=grid::gpar(fontsize=20), y= unit(0.4,"npc"))
  minus <- grid::textGrob("-",gp=grid::gpar(fontsize=20), y = unit(0.4,"npc"))
  
  library(pBrackets)
  # pBrackets dessine directement, il faut donc le transformer
  # en grob pour pouvoir le mettre dans un ggplot
  bracketsGrob <- function(...){
    l <- list(...)
    e <- new.env()
    e$l <- l
    grid:::recordGrob(  {
      do.call(pBrackets::grid.brackets, l)
    }, e)
  }
  
  accolade <- bracketsGrob(0.02,0,0.98,0,ticks=0.17,h=0.4)
  
  final <- wrap_plots(
    A,
    LA, wrap_plots(B,wrap_plots(LB,plot_spacer(),ncol=1,heights=c(0.9,0.1)),ncol=2,widths=c(3,1)),
    accolade,
    wrap_plots(C,plus,D,minus,E,ncol=5,widths = c(1,0.05,1,0.05,1)),
    design="
AAA
BCC
DDD
EEE
",heights=c(2.5,1.5,0.2,1),widths=c(2,1))
}

Fig4=globalplot(aa45_2055,MeanMeanExpositionTotPlotNorm45_2055,"A Vulnerability to global change")

pp45_2055=globalplot(aa45_2055,MeanMeanExpositionTotPlotNorm45_2055,"A Vulnerability to global change (RCP 4.5 2055)")
pp85_2055=globalplot(aa85_2055,MeanMeanExpositionTotPlotNorm85_2055,"A Vulnerability to global change (RCP 8.5 2055)")
pp45_2085=globalplot(aa45_2085,MeanMeanExpositionTotPlotNorm45_2085,"A Vulnerability to global change (RCP 4.5 2085)")
pp85_2085=globalplot(aa85_2085,MeanMeanExpositionTotPlotNorm85_2085,"A Vulnerability to global change (RCP 8.5 2085)")

ggsave("~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Rscripts/FigFutureRevision/FigFinalVulnerabilityRCP45_2055Updated.pdf",Fig4,width=5,height=6.5)


ggsave("~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Rscripts/FigFutureRevision/FigFinalVulnerabilityRCP45_2055Updated.png",pp45_2055,width=5,height=6.5)
ggsave("~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Rscripts/FigFutureRevision/FigFinalVulnerabilityRCP85_2055Updated.png",pp85_2055,width=5,height=6.5)
ggsave("~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Rscripts/FigFutureRevision/FigFinalVulnerabilityRCP45_2085Updated.png",pp45_2085,width=5,height=6.5)
ggsave("~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Rscripts/FigFutureRevision/FigFinalVulnerabilityRCP85_2085Updated.png",pp85_2085,width=5,height=6.5)





