##This is the code for the analysis presented in the paper:
##Bonfim et al. Land-use homogenization reduces the occurrence and diversity of frugivorous birds in a tropical biodiversity hotspot
##Ecological applications


library(tidyverse)
library(Hmsc)
library(phangorn)
library(viridis)
library(ape)

# You need to provide an SXY file that is provided in "Data_Bonfim et al 2023".
setwd("C:\\Users\\ferna\\Documents\\Capitulo 3\\analysis\\teste_script")


getwd()

# READING IN SXY: study design (S) and/or covariates (X) and species data (Y) 
##Save the tab "SXY" from "Data_Bonfim et al 2023" as a csv file

SXY = read.csv("SXY.csv", header= T, sep= ",") # or use read.csv2 when values are semicolon-separated
dim(SXY)
# Modify the next three lines to split your SXY file to components that relate to
# S: study design, including units of study and their possible coordinates (named as lon and lat to indicate that they relate to the units of ID_study)
# and also the Phytophysiognomy
# X: covariates to be used as predictors
# Y: species data

S=SXY[,1:6]
S$ID_study <- as.factor(S$ID_study)
S$ID_site <- as.factor(S$ID_site)
S$Phytophysiognomy <- as.factor((S$Phytophysiognomy))
str(S)


X = SXY %>%
  select(FC_5000,agriculture_5000, Pasture_1000, altitude_750, bio12_4500, count_effort)

X$count_effort <- log(X$count_effort)
str(X)

Y=as.matrix(SXY[,74:140])
dim(Y)

# Check that the data looks as it should!
View(S)
View(X)
View(Y)
# check that community data are numeric and have finite numbers. If the script
# writes "Y looks OK", you are ok.
if (is.numeric(as.matrix(Y)) || is.logical(as.matrix(Y)) && is.finite(sum(Y, na.rm=TRUE))) {
  print("Y looks OK")
} else {
  print("Y should be numeric and have finite values")	}
# Check that the stydy design data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(S))) {
  print("S has NA values - not allowed for")
} else {
  print("S looks ok")	}
# Check that the covariate data do not have missing values (they are allowed for Y but not S, X, P or Tr)
if (any(is.na(X))) {
  print("X has NA values - not allowed for")
} else {
  print("X looks ok")	}

##For the trais, use the tab "traits" in "Data_Bonfim et al 2023". Save the file as csv
TP = read.csv("traits.csv", header= T, sep= ",")
str(TP)

rownames(TP) <- colnames(Y)
# The script below checks if the species names in TP are identical and in the same order as in Y
# If the script prints "species names in TP and SXY match", you are ok.
# If it says that they do not match, you need to modify the files so that they match 
str(Y)
  if(all(rownames(TP)==colnames(Y))) {
    print("species names in TP and SXY match")
  } else{
    print("species names in TP and SXY do not match")
  }

  # Check that the Tr data do not have missing values (they are allowed for Y but not S, X, P or Tr)
  if (any(is.na(TP))) {
    print("Tr has NA values - not allowed for")
  } else {
    print("Tr looks ok")	}
  

##To use the phylogenetic tree you must get the species names in the tab
## "traits" and put in the site birdtree.org to get the trees. Use the Hackett all species set.
##After this, follow this way

##Import the tress
tree.all <- read.nexus("output.nex")

##Get a consensus tree
maxTree<-maxCladeCred(tree.all, tree = TRUE)

P = maxTree
# When you look at P (e.g. write P and press enter),
# you should see that it is a phylogenetic tree which
# is rooted and includes branch lengths and tip labels
# The script below checks if the species names in P are identical (but not necessarily in the same order) as in Y
# If the script prints "species names in P and SXY match", you are ok.
# If it says that they do not match, you need to modify the files so that they match 
if(all(sort(P$tip.label) == sort(colnames(Y)))){
  print("species names in P and SXY match")
} else{
  print("species names in P and SXY do not match")
}

  
##Fitting the model
  
##Coordinates
xy = as.matrix(cbind(S$lon, S$lat))
str(xy)
head(xy)
  
#SETTING UP THE HMSC-MODEL

  
studyDesign = data.frame(ID_study = as.factor(S$ID_study), Phytophysiognomy = as.factor(S$Phytophysiognomy))

rownames(xy)=studyDesign[,1]

rL = HmscRandomLevel(sData=xy)
str(rL)
head(rL)

rL2 <- HmscRandomLevel(units=S$Phytophysiognomy)
str(rL2)
  
##X data
XData <- data.frame(X)
str(XData)
  
XFormula = ~ FC_5000+agriculture_5000+Pasture_1000+altitude_750+bio12_4500+count_effort
  
##Trait data
rownames(TP) <- colnames(Y)
  
TrData = data.frame(TP)
str(TrData)
  
TrFormula = ~frugivory_degree+body_mass+HWI+ForStrat_ground+ForStrat_understory+ForStrat_midhigh+ForStrat_canopy
  
  
nChains = 1
samples = 5
thin = 5
  
model <- Hmsc(Y=Y, XData = XData, XFormula=XFormula, distr="probit",
                TrData = TrData, TrFormula = TrFormula, phyloTree = maxTree,
                studyDesign=studyDesign, ranLevels=list("ID_study"=rL, "Phytophysiognomy"=rL2))
  
  
m = sampleMcmc(model, samples = samples, thin = thin,
                 transient = ceiling(0.5*samples*thin),
                 nChains = nChains, nParallel = 6)


##Checking model convergence

mpost_m = convertToCodaObject(m)
gelman.diag(mpost_m$Beta,multivariate=FALSE)$psrf


##Evaluating model fitting.

preds_m = computePredictedValues(m)
result_evaluate_m <-evaluateModelFit(hM=m, predY=preds_m)

mean(result_evaluate_m$AUC)

##Computing variance partitioning
##Figures in paper were edited in slide editor

VP<- computeVariancePartitioning (m)

plotVariancePartitioning(m, VP, cols=viridis(8), main="")

##Calculating the beta parameters

postBeta = getPostEstimate(m, parName="Beta")

plotBeta(m, post=postBeta, param="Support", plotTree=T, 
         spNamesNumbers=c(T,F), supportLevel = .95,
         colors = colorRampPalette(c("red", "white", "blue")))

##Calculating the gamma parameters

postGamma = getPostEstimate(m, parName="Beta")
plotGamma(m, post=postGamma, param="Support", supportLevel = .95, 
          trNamesNumbers=c(TRUE,TRUE), colors = colorRampPalette(c("red", "white", "blue")))


##Predicting for forest cover
Gradient_fc = constructGradient(m, focalVariable="FC_5000")

pred.fc = predict(m, XData=Gradient_fc$XDataNew, studyDesign=Gradient_fc$studyDesignNew,
                  ranLevels=Gradient_fc$rLNew, expected=TRUE)

plotGradient(m, Gradient_fc, pred = pred.fc, measure="S", index=6, las=2, cex.lab = 1.9, cex.axis = 1,
             showData = TRUE, xlabel = "",  ylabel = "", main = "", showPosteriorSupport = T,
             pointsize = 0.6, pointcol = "gray30")

##Predicting for agriculture
Gradient_agri = constructGradient(m, focalVariable="agriculture_5000")

pred.agri = predict(m, XData=Gradient_agri$XDataNew, studyDesign=Gradient_agri$studyDesignNew,
                    ranLevels=Gradient_agri$rLNew, expected=TRUE)

plotGradient(m, Gradient_agri, pred=pred.agri, measure="S", index=7, las=1, cex.lab = 1.9, cex.axis = 1,
             showData = TRUE, xlabel = "",  ylabel = "", showPosteriorSupport = T,
             pointsize = .6, pointcol = "gray30")

##Predicting for pasture
Gradient_past = constructGradient(m, focalVariable="Pasture_1000")

pred.past = predict(m, XData=Gradient_past$XDataNew, studyDesign=Gradient_past$studyDesignNew,
                    ranLevels=Gradient_past$rLNew, expected=TRUE)

plotGradient(m, Gradient_past, pred=pred.past, measure="S", index = 1, las=1, cex.lab = 1.9, cex.axis = 1,
             showData = TRUE, xlabel = "",  ylabel = "", showPosteriorSupport = T,
             pointsize = .6, pointcol = "gray30")







  