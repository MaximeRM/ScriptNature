# Required packages
libs <- c("ade4","ape","SCGLR","ggplot2","raster","maptree","ggalt","viridis","data.table","scales","phylobase","adephylo","entropart")
for(l in libs) {
  library(l, character.only = TRUE)
}

# Required functions
source("Functions/TAG.R")
source("Functions/GGfuns.R")

# Load data
abond <- readRDS("Data/CoForData.rds") # Abundance data (available upon request see http://dx.doi.org/10.18167/DVN1/UCNCA7)
tabTrait <- readRDS("Data/traitTab.rds") # Trait data (available at http://dx.doi.org/10.18167/DVN1/UCNCA7)
GridTot <- readRDS("Data/GridTot.rds") # Study area with predictors (available at http://dx.doi.org/10.18167/DVN1/UCNCA7)
Anthr <- stack("Data/Anthropogenicpressure.tif")  # Human forest-disturbance index in the present and in 2085 (available at http://dx.doi.org/10.18167/DVN1/UCNCA7)

#################################################
# ORGANIZE DATA 

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

# Summary statistics for taxonomic assemblages
dudiCA <- dudi.coa(abond[,nY], scann =FALSE, nf = 10) # Correspondence analysis (CA) done on observations
compowholeCA <- suprow(dudiCA, mat_pred_val[, nY]) # Projecting predictions in the CA performed on observations 
diag(cor(dudiCA$li, compowholeCA$lisup,method="spearman"))

# Summary statistics for functional assemblages
cwmTraitobs  <-  as.data.frame(TAG(tabTrait[nY,],abond[,nY])$tag)
cwmTraitpred  <- as.data.frame(TAG(tabTrait[nY,],mat_pred_val[,nY])$tag)
diag(cor(cwmTraitobs, cwmTraitpred, method="spearman"))

#################################################
# Predicting abundances over the study area using the fitted SCGLR model

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
saveRDS(pred,file="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/Round3/Soumission/Revision/Submission/Dataverse/Predicted_abundances/PredictedAbundance.rds")



#################################################
# VULNERABILITY ANALYSIS

# Load present climate variables
dataClim=stack("Data/OUTscenarios/dataClimPresent.grd")
names(dataClim)=c(paste0("C",1:19),c("meanET0","meanCWB","sumCWD","maxCWD","MCWD"))
dataClimDat=as.data.frame(dataClim)
coord=coordinates(dataClim)
# Compute present climate components
meanPres <- genus.scglr$centerx
sdPres <- genus.scglr$invsqrtm
scaledVar <- as.matrix(sweep(dataClimDat[,nX],2,meanPres,FUN="-"))%*%sdPres
CompAfricaPres=as.data.frame(as.matrix(scaledVar)%*%as.matrix(genus.scglr$u))
names(CompAfricaPres)=c("CC1","CC2","CC3")
CompAfricaPresRast=rasterFromXYZ(cbind(coord,CompAfricaPres))

# #####################
# SENSITIVITY CLIMATE

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
TAGres <- as.data.frame(TAG(datSensitivitySPobs,prediction[,nY])$tag)
filt <- GridTot$CCI != 160
MeanSensitivityClim  <- rasterFromXYZ(cbind(prediction[filt,c("lon","lat")],1/TAGres$meanCCsd[filt]))
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
  CompFutur=as.data.frame(as.matrix(scaledVar)%*%as.matrix(genus.scglr$u))
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
listClimateModel=read.csv("Data/OUTscenarios/listClimateModel.csv")
  
resRCP45_2085=stack()
for (i in 1:nrow(listClimateModel)){
  fileFut45_85=paste0("Data/OUTscenarios/dataClimFut_",listClimateModel[i,2],"_",listClimateModel[i,1],"_RCP45_2085.grd")
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

# List of Genera in dataset
Gendata <- sapply(strsplit(nY,split="_"),"[",1)
  
# Phylogenetic tree: Janssens, S. B. et al (2020). A large-scale species level dated angiosperm phylogeny for evolutionary and ecological analyses. Biodiversity data journal, 8.
PhyloTree <- read.nexus("Data/bdj-08-e39677-s005.tre")
taxa <- PhyloTree$tip.label
taxaGen <- sapply(strsplit(taxa,split="_"),"[",1)

# Prune the tree with only the genera of the dataset
filtsp <- taxa%in%nY
TipToKeepSp <- taxa[filtsp][!duplicated(taxaGen[filtsp])]
filtgen <- taxaGen%in%Gendata & !filtsp & !taxaGen%in%taxaGen[taxa%in%nY]
TipToKeepGen <- taxa[filtgen][!duplicated(taxaGen[filtgen])]
Tree <- keep.tip(PhyloTree,tip=c(TipToKeepSp,TipToKeepGen))

# Change names to genus level
Tree$tip.label=sapply(strsplit(Tree$tip.label,split="_"),"[",1)

# Sum individuals per genus
dataPred=prediction[,nY]
dataPredt=as.data.frame(t(dataPred))
dataPredt$genus=Gendata
dataPredt=data.table(dataPredt)
dataPredtg=dataPredt[,lapply(.SD, sum, na.rm=TRUE),by=genus]
genus=dataPredtg$genus
dataPredtg$genus=NULL
dataPredg=as.data.frame(t(dataPredtg))
colnames(dataPredg)=genus
dataPred=dataPredg[,Tree$tip.label] 

# Test phylogenetic signal in component effects
# Effect of component on each taxa
taxScore=genus.scglr$gamma
taxScoret=as.data.frame(t(taxScore))
taxScoret$genus=Gendata
taxScoret=data.table(taxScoret)
taxScoretg=taxScoret[,lapply(.SD, mean, na.rm=TRUE),by=genus]
genus=taxScoretg$genus
taxScoretg$genus=NULL
taxScoretg=as.data.frame(taxScoretg)
rownames(taxScoretg)=genus
taxScore=taxScoretg[Tree$tip.label,2:4]
tree_in_phylo4 <- as(Tree,"phylo4")
aa=proxTips(tree_in_phylo4,method="oriAbouheif")
abouheif.moran(taxScore, aa,nrepet = 9999)

# Chao index with q=1
dataPredChao=t(apply(dataPred,1,function(x) as.ProbaVector(x)))
# ChaoDiv=apply(dataPredChao,1,function(x) ChaoPD(x,q=1,Tree)) # Time consuming
# saveRDS(ChaoDiv,"Data/ChaoDiv.rds")
ChaoDiv=readRDS("Data/ChaoDiv.rds")
AdaptationDat<-cbind(prediction[,c("lon","lat")],ChaoDiv)
AdaptationRast=rasterFromXYZ(AdaptationDat)
AdaptationRast[is.na(MeanSensitivityClim)]=NA
plot(AdaptationRast)

# #############################################################################
## CLIMATE VULNERABILITY 

########################################################
# Normalize the three components of vulnerability
MeanSensitivityNorm=calc(MeanSensitivityClim,scales::rescale)
MeanExpositionTotNorm45_2085=calc(MeanExpositionClim$layer,scales::rescale)
AdaptationNorm=calc(AdaptationRast,scales::rescale)

# Calculate vulnerability
VulnerabilityClim45_2085=MeanSensitivityNorm+MeanExpositionTotNorm45_2085-AdaptationNorm

VulnerabilityMaps=stack(MeanSensitivityClim,
                        MeanExpositionClim,
                        AdaptationRast,
                        VulnerabilityClim45_2085)
names(VulnerabilityMaps)=c("Sensitivity","Exposure","Adaptation","Climate_vulnerabilty")
plot(VulnerabilityMaps)

writeRaster(VulnerabilityMaps,file="~/Rejou/Post-doc_Cofortip/AnalyseCOFORTIP/SCGLR/NewMs/Nature/Submitted/Round2/Round3/Soumission/Revision/Submission/Dataverse/Climate_vulnerability/ClimateVulnerability.tif",
            format="GTiff",overwrite=TRUE)

# #############################################################################
## ANTHROPIC VULNERABILITY 

VulnerabilityAnthr=resample(Anthr$Anthropogenicpressure.2,MeanSensitivityClim)
VulnerabilityAnthr=crop(VulnerabilityAnthr,MeanSensitivityClim)
VulnerabilityAnthr[is.na(MeanSensitivityClim)]=NA
plot(VulnerabilityAnthr)

# #############################################################################
## VULNERABILITY TO GLOBAL CHANGES
library(colorplaner)# remotes::install_github("wmurphyrd/colorplaner") # the 23/04/2020

VulnerabilityAnthrNorm=calc(VulnerabilityAnthr,scales::rescale)
VulnerabilityClimNorm45_2085=calc(VulnerabilityClim45_2085,scales::rescale)
globVulnerability=stack(VulnerabilityClimNorm45_2085,VulnerabilityAnthrNorm)

datVuln=as.data.frame(globVulnerability,xy=T)
names(datVuln)=c("x","y","VulnerabilityClim","VulnerabilityAnthr")
# Apply a strech to variables
Strech <- function(x){
  v <- quantile(x, c(0.01,0.99), na.rm = TRUE)
  x[!is.na(x) & x < v[1]] <- v[1]
  x[!is.na(x) & x > v[2]] <- v[2]
  out=normalize(x)
}
datVuln$VulnerabilityClim=Strech(datVuln$VulnerabilityClim)
datVuln$VulnerabilityAnthr=Strech(datVuln$VulnerabilityAnthr)

VulnTot <- ggplot()+
  labs(x="long", y="lat")+
  coord_fixed(xlim=c(8.1,27), ylim=c(-6,5))+
  geom_tile(aes(x = x, y = y, fill = VulnerabilityClim, fill2 = VulnerabilityAnthr),data=datVuln)+
  scale_fill_colorplane(name="toto",axis_title ="Vulnerability to \nclimate change", 
                        axis_title_y = "B Anthropogenic \npressure in 2085", 
                        color_projection = YUV_projection, Y=0.65,
                        na.color = NA,
                        breaks =c(0.1,0.9), breaks_y =c(0.1,0.9), 
                        labels =c("Low","High"), labels_y =c("Low","High"))

