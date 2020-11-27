## This is the code for the manuscript "Landscape composition and configuration
##drive divergent patterns for taxonomic and functional diversity of tropical
##frugivorous birds. By Fernando César Gonçalves Bonfim, Pavel Dodonov and Eliana Cazetta##
##E-mail: fernandouesb@gmail.com##

##Load your library

setwd("C:\\Users....")

##Calculating functional diversity index
library(FD)
library(ape)
library(geometry)
library(misc3d)
library(ptinpoly)
library(vegan)
library(ade4)
library(spatstat)
library(deldir)
library(picante)
library(sp)
library(gstat)
library(MuMIn)
library(bbmle)
library(PerformanceAnalytics)
library(MASS)
library(lme4)##Adjust glmm negative binomial
library(nlme)
library(scales)
library(mctest)



##FIRST PART##

##The first part is calculate functional diversity. 
## The values calculated are in table "Analysis" in the apendix##
##If you do not want to check them, go to the second part.

##Loading files
##traits matrix ("traits" in Appendix)
traits_raw <- read.csv(file="traits.csv", header=T, row.names=1, sep="," )
str(traits_raw)

##transforming variable
traits_raw$Frugivory_score <- as.ordered(traits_raw$Frugivory_score)

##loading matrix of presence in each site ("presence_ausence" in Appendix)
presence <- read.csv(file="presence_ausence.csv", header=T, row.names=1, sep=",") 
str(presence)
presence <- as.matrix(presence)

# Distance matrix required as traits are not all continous and quantitative:

trait <- gowdis(traits_raw, ord="podani")
trait=as.matrix(trait)
str(trait)

###########################################
resOrig <- dbFD(
  trait, presence, 
  corr="none", calc.FRic=TRUE, 
  stand.FRic=FALSE, scale.RaoQ=FALSE, calc.FGR=FALSE,
  calc.CWM=FALSE, print.pco=FALSE
)
names(resOrig)

###  Null model permutations ######
set.seed(1)
null.model.type <- "independentswap" # available methods: "independentswap", "frequency", "richness", "trialswap"
nperm <- 1000
FRic <- matrix(NaN, nrow=nperm, ncol=length(resOrig$FRic))
FDis <- matrix(NaN, nrow=nperm, ncol=length(resOrig$FDis))

presence.null <- vector(mode="list", nperm)
t1 <- Sys.time()
for(i in seq(nperm)){ # loop for each permutation
  presence.i <- randomizeMatrix(presence, null.model = null.model.type)
  traits.i <- trait
  
  presence.null[[i]] <- presence.i
  
  # Record results
  res.i <- dbFD(
    traits.i, presence.i, 
    corr="none", calc.FRic=TRUE,
    stand.FRic=FALSE, scale.RaoQ=FALSE, calc.FGR=FALSE,
    calc.CWM=FALSE, print.pco=FALSE
  )
  # Save result
  FRic[i,] <- res.i$FRic
  FDis[i,] <- res.i$FDis
  rm(presence.i); rm(res.i)
  print(paste0(round(i/nperm*100), "% of permutations completed"))
}
t2 <- Sys.time()
t2 - t1 # time elapsed
(t2 - t1)/nperm * 999 # time required for 999 permutations

### Summary stats
# FRic
SESFRic = (resOrig$FRic - apply(FRic, 2, mean)) / apply(FRic, 2, sd) # standardized effect size
qFRic <- NaN*resOrig$FRic
for(i in seq(qFRic)){
  qFRic[i] <- sum(resOrig$FRic[i] > FRic[,i]) / length(FRic[,i])
}
sigFRic <- qFRic < 0.05 | qFRic > 0.95 # test if outside distribution


##SECOND PART##
##Generate the models

setwd("...")

dados<- read.table("analisys.csv", header = T, sep = ",")

str(dados)
dados$year <- as.factor(dados$year)
dados$coords<-as.matrix(cbind(dados$lon,dados$lat))

##Let's start with a full model
mix.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(mnnd)+scale(np)+coords, offset(log(dados$count_effort)),
                data = dados, na.action = "na.fail")

##Checking VIF
mc.plot(mix.all)
##Checking resisuals
par(mfrow=c(2,2))
plot(mix.all, pch=19)

##Let's construct semivariogram
E <- resid(mix.all)
mydata <- data.frame(E, dados$lon, dados$lat)
coordinates(mydata) <- c("dados.lon", "dados.lat")

Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Frugivorous richness", cex=1.5, col="black")


##Lets try binomial negative
mod.fc_ed_mnnd.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(mnnd)+coords, offset(log(dados$count_effort)),
                           data = dados, na.action = "na.fail")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords, offset(log(dados$count_effort)),
                                                 data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords, offset(log(dados$count_effort)),
                                     data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords, offset(log(dados$count_effort)),
                                       data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords, offset(log(dados$count_effort)),
                                       data = dados,  na.action = "na.fail")

mod.fc_ed.all<-glm.nb(fr~scale(fc)+scale(ed)+ coords, offset(log(dados$count_effort)),
                      data = dados, na.action = "na.fail")

mod.fc_ed_int.fc_ed.all<-glm.nb(fr~scale(fc)+scale(ed)+scale(fc):scale(ed) + coords, offset(log(dados$count_effort)),
                                data = dados, na.action = "na.fail")

mod.fc_mnnd.all<-glm.nb(fr~scale(fc)+scale(mnnd)+ coords, offset(log(dados$count_effort)),
                        data = dados, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.all<-glm.nb(fr~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords, offset(log(dados$count_effort)),
                                     data = dados, na.action = "na.fail")

mod.ed_mnnd.all<-glm.nb(fr~scale(ed)+scale(mnnd)+ coords, offset(log(dados$count_effort)),
                        data = dados, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.all<-glm.nb(fr~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords, offset(log(dados$count_effort)),
                                    data = dados, na.action = "na.fail")

mod.fc.all<-glm.nb(fr~scale(fc)+ coords, offset(log(dados$count_effort)),
                   data = dados,na.action = "na.fail")##Model do not converged

mod.ed.all<-glm.nb(fr~scale(ed)+ coords, offset(log(dados$count_effort)),
                   data = dados, na.action = "na.fail")

mod.mnnd.all<-glm.nb(fr~scale(mnnd)+ coords, offset(log(dados$count_effort)),
                     data = dados, na.action = "na.fail")

##To calculate variable importance

wAICc.all <- AICctab(mod.fc_ed_mnnd.all, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.all, mod.fc_ed_mnnd_int.fc_ed.all,
                     mod.fc_ed_mnnd_int.fc_mnnd.all, mod.fc_ed_mnnd_int.mnnd_ed.all, mod.fc_ed.all,
                     mod.fc_ed_int.fc_ed.all, mod.fc_mnnd.all, mod.fc_mnnd_int.fc_mnnnd.all,
                     mod.ed_mnnd.all,mod.ed_mnnd_int.ed_mnnd.all, mod.fc.all, mod.ed.all, mod.mnnd.all, weights=T)

wAICc.all
summary(mod.fc_ed_int.fc_ed.all)

nomes.modelos.all <- attr(wAICc.all, "row.names")

wAICc.modelos.all <- wAICc.all$weight

modelos.fc.all <- grep("fc", nomes.modelos.all)

modelos.ed.all <- grep("ed", nomes.modelos.all)

modelos.mnnd.all <- grep("mnnd", nomes.modelos.all)

importance.fc.all <- sum(wAICc.modelos.all[modelos.fc.all])
importance.ed.all <- sum(wAICc.modelos.all[modelos.ed.all])
importance.mnnd.all <- sum(wAICc.modelos.all[modelos.mnnd.all])

importancias.all <- c(importance.fc.all, importance.ed.all, importance.mnnd.all)
names(importancias.all) <- c("fc", "ed", "mnnd")
importancias.all


##To calculate mean effect size
model.all<-(model.avg(mod.fc_ed_mnnd.all, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.all, mod.fc_ed_mnnd_int.fc_ed.all,
                      mod.fc_ed_mnnd_int.fc_mnnd.all, mod.fc_ed_mnnd_int.mnnd_ed.all, mod.fc_ed.all,
                      mod.fc_ed_int.fc_ed.all, mod.fc_mnnd.all, mod.fc_mnnd_int.fc_mnnnd.all,
                      mod.ed_mnnd.all, mod.ed_mnnd_int.ed_mnnd.all, mod.fc.all, mod.ed.all, mod.mnnd.all))
summary(model.all)

####Fric
##Because for functional diversity we excluded 12 sites, we need to exclude lines where there are NA's in SESFRic and FDis
dados$SESFRic
dados_fun<-dados2<-dados[-c(18, 19, 21, 23, 79, 83, 94, 109, 126, 127, 130, 131),]
dados_fun$year <- as.factor(dados_fun$year)
str(dados_fun)
dados_fun$coords<-as.matrix(cbind(dados_fun$lon,dados_fun$lat))


mod.fc_ed_mnnd.fric<-lm(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ coords,
                        data = dados_fun, na.action = "na.fail")

par(mfrow=c(2,2))
plot(mod.fc_ed_mnnd.fric, pch=19)



E <- resid(mod.fc_ed_mnnd.fric)
mydata <- data.frame(E, dados_fun$lon, dados_fun$lat)
coordinates(mydata) <- c("dados_fun.lon", "dados_fun.lat")

Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "SES Functional richness", cex=1.5, col="black")


###Lets construct the models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fric<-lm(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords,
                                              data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.fric<-lm(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords,
                                  data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.fric<-lm(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords,
                                    data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.fric<-lm(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                    data = dados_fun, na.action = "na.fail")

mod.fc_ed.fric<-lm(SESFRic~scale(fc)+scale(ed)+ coords,
                   data = dados_fun, na.action = "na.fail")

mod.fc_ed_int.fc_ed.fric<-lm(SESFRic~scale(fc)+scale(ed)+scale(fc):scale(ed) + coords,
                             data = dados_fun, na.action = "na.fail")

mod.fc_mnnd.fric<-lm(SESFRic~scale(fc)+scale(mnnd)+ coords,
                     data = dados_fun, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.fric<-lm(SESFRic~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords,
                                  data = dados_fun, na.action = "na.fail")

mod.ed_mnnd.fric<-lm(SESFRic~scale(ed)+scale(mnnd)+ coords,
                     data = dados_fun, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.fric<-lm(SESFRic~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                 data = dados_fun, na.action = "na.fail")

mod.fc.fric<-lm(SESFRic~scale(fc)+ coords,
                data = dados_fun, na.action = "na.fail")

mod.ed.fric<-lm(SESFRic~scale(ed)+ coords,
                data = dados_fun, na.action = "na.fail")

mod.mnnd.fric<-lm(SESFRic~scale(mnnd)+ coords,
                  data = dados_fun, na.action = "na.fail")


##To calculate variable importance

wAICc.fric <- AICctab(mod.fc_ed_mnnd.fric, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fric, mod.fc_ed_mnnd_int.fc_ed.fric,
                      mod.fc_ed_mnnd_int.fc_mnnd.fric, mod.fc_ed_mnnd_int.mnnd_ed.fric, mod.fc_ed.fric,
                      mod.fc_ed_int.fc_ed.fric, mod.fc_mnnd.fric, mod.fc_mnnd_int.fc_mnnnd.fric,
                      mod.ed_mnnd.fric,mod.ed_mnnd_int.ed_mnnd.fric, mod.fc.fric, mod.ed.fric, mod.mnnd.fric, weights=T)

wAICc.fric
summary(mod.ed.fric)

nomes.modelos.fric <- attr(wAICc.fric, "row.names")

wAICc.modelos.fric <- wAICc.fric$weight

modelos.fc.fric <- grep("fc", nomes.modelos.fric)

modelos.ed.fric <- grep("ed", nomes.modelos.fric)

modelos.mnnd.fric <- grep("mnnd", nomes.modelos.fric)

importance.fc.fric <- sum(wAICc.modelos.fric[modelos.fc.fric])
importance.ed.fric <- sum(wAICc.modelos.fric[modelos.ed.fric])
importance.mnnd.fric <- sum(wAICc.modelos.fric[modelos.mnnd.fric])

importancias.fric <- c(importance.fc.fric, importance.ed.fric, importance.mnnd.fric)
names(importancias.fric) <- c("fc", "ed", "mnnd")
importancias.fric

##To calculate mean effect size
model.fric<-(model.avg(mod.fc_ed_mnnd.fric, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fric, mod.fc_ed_mnnd_int.fc_ed.fric,
                       mod.fc_ed_mnnd_int.fc_mnnd.fric, mod.fc_ed_mnnd_int.mnnd_ed.fric, mod.fc_ed.fric,
                       mod.fc_ed_int.fc_ed.fric, mod.fc_mnnd.fric, mod.fc_mnnd_int.fc_mnnnd.fric,
                       mod.ed_mnnd.fric,mod.ed_mnnd_int.ed_mnnd.fric, mod.fc.fric, mod.ed.fric, mod.mnnd.fric))
summary(model.fric)


####FDis
mod.fc_ed_mnnd.fdis<-lm(FDis~scale(fc)+scale(ed)+scale(mnnd)+ coords,
                        data = dados_fun, na.action = "na.fail")

par(mfrow=c(2,2))
plot(mod.fc_ed_mnnd.fdis, pch=19)

E <- resid(mod.fc_ed_mnnd.fdis)
mydata <- data.frame(E, dados_fun$lon, dados_fun$lat)
coordinates(mydata) <- c("dados_fun.lon", "dados_fun.lat")

Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "FDis", cex=1.5, col="black")


###Lets construct the models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fdis<-lm(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords,
                                              data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.fdis<-lm(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords,
                                  data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.fdis<-lm(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords, 
                                    data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.fdis<-lm(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                    data = dados_fun, na.action = "na.fail")

mod.fc_ed.fdis<-lm(FDis~scale(fc)+scale(ed)+ coords,
                   data = dados_fun, na.action = "na.fail")

mod.fc_ed_int.fc_ed.fdis<-lm(FDis~scale(fc)+scale(ed)+scale(fc):scale(ed) + coords,
                             data = dados_fun, na.action = "na.fail")

mod.fc_mnnd.fdis<-lm(FDis~scale(fc)+scale(mnnd)+ coords,
                     data = dados_fun, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.fdis<-lm(FDis~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords,
                                  data = dados_fun, na.action = "na.fail")

mod.ed_mnnd.fdis<-lm(FDis~scale(ed)+scale(mnnd)+ coords,
                     data = dados_fun, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.fdis<-lm(FDis~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                 data = dados_fun, na.action = "na.fail")

mod.fc.fdis<-lm(FDis~scale(fc)+ coords,
                data = dados_fun, na.action = "na.fail")

mod.ed.fdis<-lm(FDis~scale(ed)+ coords,
                data = dados_fun, na.action = "na.fail")

mod.mnnd.fdis<-lm(FDis~scale(mnnd)+ coords,
                  data = dados_fun, na.action = "na.fail")


wAICc.fdis <- AICctab(mod.fc_ed_mnnd.fdis, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fdis, mod.fc_ed_mnnd_int.fc_ed.fdis,
                      mod.fc_ed_mnnd_int.fc_mnnd.fdis, mod.fc_ed_mnnd_int.mnnd_ed.fdis, mod.fc_ed.fdis,
                      mod.fc_ed_int.fc_ed.fdis, mod.fc_mnnd.fdis, mod.fc_mnnd_int.fc_mnnnd.fdis,
                      mod.ed_mnnd.fdis,mod.ed_mnnd_int.ed_mnnd.fdis, mod.fc.fdis, mod.ed.fdis, mod.mnnd.fdis, weights=T)
wAICc.fdis


nomes.modelos.fdis <- attr(wAICc.fdis, "row.names")

wAICc.modelos.fdis<- wAICc.fdis$weight

modelos.fc.fdis <- grep("fc", nomes.modelos.fdis)

modelos.ed.fdis <- grep("ed", nomes.modelos.fdis)

modelos.mnnd.fdis <- grep("mnnd", nomes.modelos.fdis)

importance.fc.fdis <- sum(wAICc.modelos.fdis[modelos.fc.fdis])
importance.ed.fdis <- sum(wAICc.modelos.fdis[modelos.ed.fdis])
importance.mnnd.fdis <- sum(wAICc.modelos.fdis[modelos.mnnd.fdis])

importancias.fdis <- c(importance.fc.fdis, importance.ed.fdis, importance.mnnd.fdis)
names(importancias.fdis) <- c("fc", "ed", "mnnd")
importancias.fdis

##To calculate effect size

model.fdis<-(model.avg(mod.fc_ed_mnnd.fdis, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fdis, mod.fc_ed_mnnd_int.fc_ed.fdis,
                       mod.fc_ed_mnnd_int.fc_mnnd.fdis, mod.fc_ed_mnnd_int.mnnd_ed.fdis, mod.fc_ed.fdis,
                       mod.fc_ed_int.fc_ed.fdis, mod.fc_mnnd.fdis, mod.fc_mnnd_int.fc_mnnnd.fdis,
                       mod.ed_mnnd.fdis,mod.ed_mnnd_int.ed_mnnd.fdis, mod.fc.fdis, mod.ed.fdis, mod.mnnd.fdis))
summary(model.fdis)


###########################################

##Body mass
str(dados)

mod.fc_ed_mnnd.body<-lm(Body~scale(fc)+scale(ed)+scale(mnnd)+ coords,
                        data = dados, na.action = "na.fail")

par(mfrow=c(2,2))
plot(mod.fc_ed_mnnd.body, pch=19)

E <- resid(mod.fc_ed_mnnd.body)
mydata <- data.frame(E, dados$lon, dados$lat)
coordinates(mydata) <- c("dados.lon", "dados.lat")

Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Body mass", cex=1.5, col="black")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.body<-lm(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords,
                                              data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.body<-lm(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords,
                                  data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.body<-lm(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords,
                                    data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.body<-lm(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                    data = dados, na.action = "na.fail")

mod.fc_ed.body<-lm(Body~scale(fc)+scale(ed)+ coords,
                   data = dados, na.action = "na.fail")

mod.fc_ed_int.fc_ed.body<-lm(Body~scale(fc)+scale(ed)+scale(fc):scale(ed) + coords,
                             data = dados, na.action = "na.fail")

mod.fc_mnnd.body<-lm(Body~scale(fc)+scale(mnnd)+ coords,
                     data = dados, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.body<-lm(Body~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords,
                                  data = dados, na.action = "na.fail")

mod.ed_mnnd.body<-lm(Body~scale(ed)+scale(mnnd)+ coords,
                     data = dados, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.body<-lm(Body~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                 data = dados, na.action = "na.fail")

mod.fc.body<-lm(Body~scale(fc)+ coords,
                data = dados, na.action = "na.fail")

mod.ed.body<-lm(Body~scale(ed)+ coords,
                data = dados, na.action = "na.fail")

mod.mnnd.body<-lm(Body~scale(mnnd)+ coords,
                  data = dados, na.action = "na.fail")


wAICc.body <- AICctab(mod.fc_ed_mnnd.body, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.body, mod.fc_ed_mnnd_int.fc_ed.body,
                      mod.fc_ed_mnnd_int.fc_mnnd.body, mod.fc_ed_mnnd_int.mnnd_ed.body, mod.fc_ed.body,
                      mod.fc_ed_int.fc_ed.body, mod.fc_mnnd.body, mod.fc_mnnd_int.fc_mnnnd.body,
                      mod.ed_mnnd.body,mod.ed_mnnd_int.ed_mnnd.body, mod.fc.body, mod.ed.body, mod.mnnd.body, weights=T)
wAICc.body

nomes.modelos.body <- attr(wAICc.body, "row.names")

wAICc.modelos.body<- wAICc.body$weight

modelos.fc.body <- grep("fc", nomes.modelos.body)

modelos.ed.body <- grep("ed", nomes.modelos.body)

modelos.mnnd.body <- grep("mnnd", nomes.modelos.body)

importance.fc.body <- sum(wAICc.modelos.body[modelos.fc.body])
importance.ed.body <- sum(wAICc.modelos.body[modelos.ed.body])
importance.mnnd.body <- sum(wAICc.modelos.body[modelos.mnnd.body])

importancias.body <- c(importance.fc.body, importance.ed.body, importance.mnnd.body)
names(importancias.body) <- c("fc", "ed", "mnnd")
importancias.body

##To calculate effect size

model.body<-(model.avg(mod.fc_ed_mnnd.body, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.body, mod.fc_ed_mnnd_int.fc_ed.body,
                       mod.fc_ed_mnnd_int.fc_mnnd.body, mod.fc_ed_mnnd_int.mnnd_ed.body, mod.fc_ed.body,
                       mod.fc_ed_int.fc_ed.body, mod.fc_mnnd.body, mod.fc_mnnd_int.fc_mnnnd.body,
                       mod.ed_mnnd.body,mod.ed_mnnd_int.ed_mnnd.body, mod.fc.body, mod.ed.body, mod.mnnd.body))
summary(model.body)


###########################################

##Bill width
mod.fc_ed_mnnd.bill<-lm(Bill~scale(fc)+scale(ed)+scale(mnnd)+ coords,
                        data = dados, na.action = "na.fail")
par(mfrow=c(2,2))
plot(mod.fc_ed_mnnd.bill, pch=19)

E <- resid(mod.fc_ed_mnnd.bill)
mydata <- data.frame(E, dados$lon, dados$lat)
coordinates(mydata) <- c("dados.lon", "dados.lat")

Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Bill width", cex=1.5, col="black")


###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.bill<-lm(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords,
                                              data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.bill<-lm(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords,
                                  data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.bill<-lm(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords,
                                    data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.bill<-lm(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                    data = dados, na.action = "na.fail")

mod.fc_ed.bill<-lm(Bill~scale(fc)+scale(ed)+ coords,
                   data = dados, na.action = "na.fail")

mod.fc_ed_int.fc_ed.bill<-lm(Bill~scale(fc)+scale(ed)+scale(fc):scale(ed) + coords,
                             data = dados, na.action = "na.fail")

mod.fc_mnnd.bill<-lm(Bill~scale(fc)+scale(mnnd)+ coords,
                     data = dados, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.bill<-lm(Bill~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords,
                                  data = dados, na.action = "na.fail")

mod.ed_mnnd.bill<-lm(Bill~scale(ed)+scale(mnnd)+ coords,
                     data = dados, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.bill<-lm(Bill~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                 data = dados, na.action = "na.fail")

mod.fc.bill<-lm(Bill~scale(fc)+ coords,
                data = dados, na.action = "na.fail")

mod.ed.bill<-lm(Bill~scale(ed)+ coords,
                data = dados, na.action = "na.fail")

mod.mnnd.bill<-lm(Bill~scale(mnnd)+ coords,
                  data = dados, na.action = "na.fail")


wAICc.bill <- AICctab(mod.fc_ed_mnnd.bill, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.bill, mod.fc_ed_mnnd_int.fc_ed.bill,
                      mod.fc_ed_mnnd_int.fc_mnnd.bill, mod.fc_ed_mnnd_int.mnnd_ed.bill, mod.fc_ed.bill,
                      mod.fc_ed_int.fc_ed.bill, mod.fc_mnnd.bill, mod.fc_mnnd_int.fc_mnnnd.bill,
                      mod.ed_mnnd.bill,mod.ed_mnnd_int.ed_mnnd.bill, mod.fc.bill, mod.ed.bill, mod.mnnd.bill, weights=T)
wAICc.bill


nomes.modelos.bill <- attr(wAICc.bill, "row.names")

wAICc.modelos.bill<- wAICc.bill$weight

modelos.fc.bill <- grep("fc", nomes.modelos.bill)

modelos.ed.bill <- grep("ed", nomes.modelos.bill)

modelos.mnnd.bill <- grep("mnnd", nomes.modelos.bill)

importance.fc.bill <- sum(wAICc.modelos.bill[modelos.fc.bill])
importance.ed.bill <- sum(wAICc.modelos.bill[modelos.ed.bill])
importance.mnnd.bill <- sum(wAICc.modelos.bill[modelos.mnnd.bill])

importancias.bill <- c(importance.fc.bill, importance.ed.bill, importance.mnnd.bill)
names(importancias.bill) <- c("fc", "ed", "mnnd")
importancias.bill

##To calculate effect size

model.bill<-(model.avg(mod.fc_ed_mnnd.bill, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.bill, mod.fc_ed_mnnd_int.fc_ed.bill,
                       mod.fc_ed_mnnd_int.fc_mnnd.bill, mod.fc_ed_mnnd_int.mnnd_ed.bill, mod.fc_ed.bill,
                       mod.fc_ed_int.fc_ed.bill, mod.fc_mnnd.bill, mod.fc_mnnd_int.fc_mnnnd.bill,
                       mod.ed_mnnd.bill,mod.ed_mnnd_int.ed_mnnd.bill, mod.fc.bill, mod.ed.bill, mod.mnnd.bill))
summary(model.bill)


###########################################
##We need to exclue one line where there is no value for hand wing index
dados2<-dados[-c(130),]


dados2$coords<-as.matrix(cbind(dados2$lon,dados2$lat))

mod.fc_ed_mnnd.hand<-lm(Hand~scale(fc)+scale(ed)+scale(mnnd)+ coords,
                        data = dados2, na.action = "na.fail")
par(mfrow=c(2,2))
plot(mod.fc_ed_mnnd.hand, pch=19)

E <- resid(mod.fc_ed_mnnd.hand)
mydata <- data.frame(E, dados2$lon, dados2$lat)
coordinates(mydata) <- c("dados2.lon", "dados2.lat")

Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Hand wing index", cex=1.5, col="black")


###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.hand<-lm(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords,
                                              data = dados2, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.hand<-lm(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords,
                                  data = dados2, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.hand<-lm(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords,
                                    data = dados2, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.hand<-lm(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                    data = dados2, na.action = "na.fail")

mod.fc_ed.hand<-lm(Hand~scale(fc)+scale(ed)+ coords,
                   data = dados2, na.action = "na.fail")

mod.fc_ed_int.fc_ed.hand<-lm(Hand~scale(fc)+scale(ed)+scale(fc):scale(ed) + coords,
                             data = dados2, na.action = "na.fail")

mod.fc_mnnd.hand<-lm(Hand~scale(fc)+scale(mnnd)+ coords,
                     data = dados2, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.hand<-lm(Hand~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords,
                                  data = dados2, na.action = "na.fail")

mod.ed_mnnd.hand<-lm(Hand~scale(ed)+scale(mnnd)+ coords,
                     data = dados2, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.hand<-lm(Hand~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords,
                                 data = dados2, na.action = "na.fail")

mod.fc.hand<-lm(Hand~scale(fc)+ coords,
                data = dados2, na.action = "na.fail")

mod.ed.hand<-lm(Hand~scale(ed)+ coords,
                data = dados2, na.action = "na.fail")

mod.mnnd.hand<-lm(Hand~scale(mnnd)+ coords,
                  data = dados2, na.action = "na.fail")


wAICc.hand <- AICctab(mod.fc_ed_mnnd.hand, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.hand, mod.fc_ed_mnnd_int.fc_ed.hand,
                      mod.fc_ed_mnnd_int.fc_mnnd.hand, mod.fc_ed_mnnd_int.mnnd_ed.hand, mod.fc_ed.hand,
                      mod.fc_ed_int.fc_ed.hand, mod.fc_mnnd.hand, mod.fc_mnnd_int.fc_mnnnd.hand,
                      mod.ed_mnnd.hand,mod.ed_mnnd_int.ed_mnnd.hand, mod.fc.hand, mod.ed.hand, mod.mnnd.hand, weights=T)
wAICc.hand

nomes.modelos.hand <- attr(wAICc.hand, "row.names")

wAICc.modelos.hand<- wAICc.hand$weight

modelos.fc.hand <- grep("fc", nomes.modelos.hand)

modelos.ed.hand <- grep("ed", nomes.modelos.hand)

modelos.mnnd.hand <- grep("mnnd", nomes.modelos.hand)

importance.fc.hand <- sum(wAICc.modelos.hand[modelos.fc.hand])
importance.ed.hand <- sum(wAICc.modelos.hand[modelos.ed.hand])
importance.mnnd.hand <- sum(wAICc.modelos.hand[modelos.mnnd.hand])

importancias.hand <- c(importance.fc.hand, importance.ed.hand, importance.mnnd.hand)
names(importancias.hand) <- c("fc", "ed", "mnnd")
importancias.hand

##To calculate effect size

model.hand<-(model.avg(mod.fc_ed_mnnd.hand, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.hand, mod.fc_ed_mnnd_int.fc_ed.hand,
                       mod.fc_ed_mnnd_int.fc_mnnd.hand, mod.fc_ed_mnnd_int.mnnd_ed.hand, mod.fc_ed.hand,
                       mod.fc_ed_int.fc_ed.hand, mod.fc_mnnd.hand, mod.fc_mnnd_int.fc_mnnnd.hand,
                       mod.ed_mnnd.hand,mod.ed_mnnd_int.ed_mnnd.hand, mod.fc.hand, mod.ed.hand, mod.mnnd.hand))
summary(model.hand)




