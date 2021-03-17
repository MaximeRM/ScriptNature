# Required packages
libs <- c("ade4","SCGLR","ggplot2","raster","maptree","ggalt","viridis")
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

################
# Formula
form <- multivariateFormula(nY, c(nX), A = nT)
# Type of family distribution
family <- rep("poisson", length(nY))

################
# Crossvalidation to select the number of climate components: Time consuming!
cv <- scglrCrossVal(form, data=abond,K=6, offset=abond$offset, nfolds=10,
                    family=family, crit=list(tol=1e-6,maxit=1000),
                    method=methodSR(l=1,s=0.1,bailout = 1000))
mean.crit <- colMeans(log(cv))
kOpt <- which.min(mean.crit)-1 # kOpt=3 in the present study

################
# Run scglr on the whole dataset
genus.scglr <-scglr(form, data = abond,K = kOpt, offset = abond$offset,
                    family = family, crit = list(tol = 1e-6, maxit = 1000),
                    method = methodSR(l = 1,s = 0.1,bailout = 1000,maxiter = 100))
plot(genus.scglr,threshold=0.1,labels.size=0.6,labels.auto=FALSE,plane=c(1,2))

################
# Assess prediction error with spatial cross validations: Time consuming!

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
prediction<- data.frame(prediction, lon = GridTot$lon,lat = GridTot$lat)

# Floristic analysis 
datFlo <- prediction[,nY]
dudiPred <- dudi.coa(datFlo, scannf=FALSE, nf=5)

# Map first three predicted CA axes (Fig 1 of the paper)
dataCA <- data.frame(prediction[, c("lon","lat")],dudiPred$li[,1:3])
rastRGBflo <- rasterFromXYZ(dataCA)
plot(rastRGBflo)

# Map functional traits (TO be done)


# Build floristic groups
Naxes <- 5
TreeClust <- hclust(dist(dudiPred$li[,1:Naxes]),method = "ward.D")
SubTree <- clip.clust(TreeClust,k = ngroup,data=prediction)
extract <- function(cl, m){
  res <- data.frame(prediction[as.integer(m), c("lon","lat")], clFac=cl)
}
classPredict <- mapply(extract, names(SubTree$membership), SubTree$membership, SIMPLIFY = FALSE)
classPredict <- do.call(rbind, classPredict)

rasttest <- rasterFromXYZ(classPredict)
plot(rasttest)

