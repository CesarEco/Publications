## This is the code for the manuscript "Landscape composition taxonomic and functional diversity of tropical
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

##import the data (table "analyses" in Appendix)

dados<- read.table("dados_novos.csv", header = T, sep = ",")

str(dados)
dados$ID_study <- as.factor(dados$ID_study)
dados$coords<-as.matrix(cbind(dados$lon,dados$lat))
hist(dados$fc)
count_effort<- dados$count_effort
mod.fc_ed_mnnd.all<-glmer(fr~scale(fc)+scale(ed)+scale(mnnd)+ coords+ (1|ID_study)+ offset(log(count_effort)),
                          data = dados, family = poisson, na.action = "na.fail")


####Calculating semivariogram
E <- resid(mod.fc_ed_mnnd.all)
mydata <- data.frame(E, dados$lon, dados$lat)
coordinates(mydata) <- c("dados.lon", "dados.lat")

###semivariograma####
Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Frugivorous richness", cex=1.5, col="black")

###Lets construct models

mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.all<-glmer(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                                                data = dados,family = poisson, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.all<-glmer(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords+(1|ID_study)+ offset(log(dados$count_effort)),
                                    data = dados,family = poisson, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.all<-glmer(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords+(1|ID_study)+ offset(log(dados$count_effort)),
                                      data = dados, family= poisson, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.all<-glmer(fr~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                                      data = dados, family=poisson, na.action = "na.fail")

mod.fc_ed.all<-glmer(fr~scale(fc)+scale(ed)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                     data = dados, family=poisson, na.action = "na.fail")

mod.fc_ed_int.fc_ed.all<-glmer(fr~scale(fc)+scale(ed)+scale(fc):scale(ed)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                               data = dados, family = poisson, na.action = "na.fail")

mod.fc_mnnd.all<-glmer(fr~scale(fc)+scale(mnnd)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                       data = dados, family= poisson, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.all<-glmer(fr~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                                    data = dados, family= poisson, na.action = "na.fail")

mod.ed_mnnd.all<-glmer(fr~scale(ed)+scale(mnnd)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                       data = dados, family= poisson, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.all<-glmer(fr~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                                   data = dados, family= poisson, na.action = "na.fail")

mod.fc.all<-glmer(fr~scale(fc)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                  data = dados, family= poisson, na.action = "na.fail")

mod.ed.all<-glmer(fr~scale(ed)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                  data = dados, family= poisson, na.action = "na.fail")

mod.mnnd.all<-glmer(fr~scale(mnnd)+ coords+ (1|ID_study)+ offset(log(dados$count_effort)),
                    data = dados, family=poisson, na.action = "na.fail")

##To calculate variable importance

wAICc.all <- AICctab(mod.fc_ed_mnnd.all, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.all, mod.fc_ed_mnnd_int.fc_ed.all,
                     mod.fc_ed_mnnd_int.fc_mnnd.all, mod.fc_ed_mnnd_int.mnnd_ed.all, mod.fc_ed.all,
                     mod.fc_ed_int.fc_ed.all, mod.fc_mnnd.all, mod.fc_mnnd_int.fc_mnnnd.all,
                     mod.ed_mnnd.all,mod.ed_mnnd_int.ed_mnnd.all, mod.fc.all, mod.ed.all, mod.mnnd.all, weights=T)

wAICc.all

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


##For Fric and FDis, exclude the lines that has no data
########################################################
####Fric
dados_fun<- read.table("dados_novos_fun.csv", header = T, sep = ",")
dados_fun$ID_study <- as.factor(dados_fun$ID_study)
str(dados_fun)
dados_fun$coords<-as.matrix(cbind(dados_fun$lon,dados_fun$lat))


mod.fc_ed_mnnd.fric<-lmer(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ coords+(1|ID_study),
                          data = dados_fun, na.action = "na.fail")

####Calculating semivariogram

E <- resid(mod.fc_ed_mnnd.fric)
mydata <- data.frame(E, dados_fun$lon, dados_fun$lat)
coordinates(mydata) <- c("dados_fun.lon", "dados_fun.lat")

###semivariograma####
Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "SES Functional richness", cex=1.5, col="black")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fric<-lmer(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords+(1|ID_study),
                                                data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.fric<-lmer(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords+(1|ID_study),
                                    data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.fric<-lmer(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords+(1|ID_study),
                                      data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.fric<-lmer(SESFRic~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+(1|ID_study),
                                      data = dados_fun, na.action = "na.fail")

mod.fc_ed.fric<-lmer(SESFRic~scale(fc)+scale(ed) +coords+(1|ID_study),
                     data = dados_fun, na.action = "na.fail")

mod.fc_ed_int.fc_ed.fric<-lmer(SESFRic~scale(fc)+scale(ed)+scale(fc):scale(ed)+ coords+(1|ID_study),
                               data = dados_fun, na.action = "na.fail")

mod.fc_mnnd.fric<-lmer(SESFRic~scale(fc)+scale(mnnd)+ coords+(1|ID_study),
                       data = dados_fun, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.fric<-lmer(SESFRic~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords+(1|ID_study),
                                    data = dados_fun, na.action = "na.fail")

mod.ed_mnnd.fric<-lmer(SESFRic~scale(ed)+scale(mnnd)+ coords+(1|ID_study),
                       data = dados_fun, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.fric<-lmer(SESFRic~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+(1|ID_study),
                                   data = dados_fun, na.action = "na.fail")

mod.fc.fric<-lmer(SESFRic~scale(fc)+ coords+(1|ID_study),
                  data = dados_fun, na.action = "na.fail")

mod.ed.fric<-lmer(SESFRic~scale(ed)+ coords+(1|ID_study),
                  data = dados_fun, na.action = "na.fail")

mod.mnnd.fric<-lmer(SESFRic~scale(mnnd)+ coords+(1|ID_study),
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

#######################################################3
####FDis
mod.fc_ed_mnnd.fdis<-lmer(FDis~scale(fc)+scale(ed)+scale(mnnd)+ coords+(1|ID_study),
                          data = dados_fun, na.action = "na.fail")


####Calculating semivariogram
E <- resid(mod.fc_ed_mnnd.fdis)
mydata <- data.frame(E, dados_fun$lon, dados_fun$lat)
coordinates(mydata) <- c("dados_fun.lon", "dados_fun.lat")

###semivariograma####
Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "FDis", cex=1.5, col="black")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fdis<-lmer(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords+(1|ID_study),
                                                data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.fdis<-lmer(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords+(1|ID_study),
                                    data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.fdis<-lmer(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords+(1|ID_study), 
                                      data = dados_fun, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.fdis<-lmer(FDis~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+(1|ID_study),
                                      data = dados_fun, na.action = "na.fail")

mod.fc_ed.fdis<-lmer(FDis~scale(fc)+scale(ed)+ coords+(1|ID_study),
                     data = dados_fun, na.action = "na.fail")

mod.fc_ed_int.fc_ed.fdis<-lmer(FDis~scale(fc)+scale(ed)+scale(fc):scale(ed)+ coords+(1|ID_study),
                               data = dados_fun, na.action = "na.fail")

mod.fc_mnnd.fdis<-lmer(FDis~scale(fc)+scale(mnnd)+ coords+(1|ID_study),
                       data = dados_fun, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.fdis<-lmer(FDis~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords+(1|ID_study),
                                    data = dados_fun, na.action = "na.fail")

mod.ed_mnnd.fdis<-lmer(FDis~scale(ed)+scale(mnnd)+ coords+(1|ID_study),
                       data = dados_fun, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.fdis<-lmer(FDis~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+(1|ID_study),
                                   data = dados_fun, na.action = "na.fail")

mod.fc.fdis<-lmer(FDis~scale(fc)+ coords+(1|ID_study),
                  data = dados_fun, na.action = "na.fail")

mod.ed.fdis<-lmer(FDis~scale(ed)+ coords+(1|ID_study),
                  data = dados_fun, na.action = "na.fail")

mod.mnnd.fdis<-lmer(FDis~scale(mnnd)+ coords+(1|ID_study),
                    data = dados_fun, na.action = "na.fail")


wAICc.fdis <- AICctab(mod.fc_ed_mnnd.fdis, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.fdis, mod.fc_ed_mnnd_int.fc_ed.fdis,
                      mod.fc_ed_mnnd_int.fc_mnnd.fdis, mod.fc_ed_mnnd_int.mnnd_ed.fdis, mod.fc_ed.fdis,
                      mod.fc_ed_int.fc_ed.fdis, mod.fc_mnnd.fdis, mod.fc_mnnd_int.fc_mnnnd.fdis,
                      mod.ed_mnnd.fdis,mod.ed_mnnd_int.ed_mnnd.fdis, mod.fc.fdis, mod.ed.fdis, mod.mnnd.fdis, weights=T)
wAICc.fdis

summary(mod.ed.fdis)

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

mod.fc_ed_mnnd.body<-lmer(Body~scale(fc)+scale(ed)+scale(mnnd)+coords+ (1|ID_study),
                          data = dados, na.action = "na.fail")

####Calculating semivariogram

E <- resid(mod.fc_ed_mnnd.body)
mydata <- data.frame(E, dados$lon, dados$lat)
coordinates(mydata) <- c("dados.lon", "dados.lat")

###semivariograma####
Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Body mass", cex=1.5, col="black")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.body<-lmer(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                                data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.body<-lmer(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords+ (1|ID_study),
                                    data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.body<-lmer(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                      data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.body<-lmer(Body~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study),
                                      data = dados, na.action = "na.fail")

mod.fc_ed.body<-lmer(Body~scale(fc)+scale(ed)+ coords+ (1|ID_study),
                     data = dados, na.action = "na.fail")

mod.fc_ed_int.fc_ed.body<-lmer(Body~scale(fc)+scale(ed)+scale(fc):scale(ed)+ coords+ (1|ID_study),
                               data = dados, na.action = "na.fail")

mod.fc_mnnd.body<-lmer(Body~scale(fc)+scale(mnnd)+ coords+ (1|ID_study),
                       data = dados, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.body<-lmer(Body~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                    data = dados, na.action = "na.fail")

mod.ed_mnnd.body<-lmer(Body~scale(ed)+scale(mnnd)+ coords+ (1|ID_study),
                       data = dados, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.body<-lmer(Body~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study),
                                   data = dados, na.action = "na.fail")

mod.fc.body<-lmer(Body~scale(fc)+ coords+ (1|ID_study),
                  data = dados, na.action = "na.fail")

mod.ed.body<-lmer(Body~scale(ed)+ coords+ (1|ID_study),
                  data = dados, na.action = "na.fail")

mod.mnnd.body<-lmer(Body~scale(mnnd)+ coords+ (1|ID_study),
                    data = dados, na.action = "na.fail")


wAICc.body <- AICctab(mod.fc_ed_mnnd.body, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.body, mod.fc_ed_mnnd_int.fc_ed.body,
                      mod.fc_ed_mnnd_int.fc_mnnd.body, mod.fc_ed_mnnd_int.mnnd_ed.body, mod.fc_ed.body,
                      mod.fc_ed_int.fc_ed.body, mod.fc_mnnd.body, mod.fc_mnnd_int.fc_mnnnd.body,
                      mod.ed_mnnd.body,mod.ed_mnnd_int.ed_mnnd.body, mod.fc.body, mod.ed.body, mod.mnnd.body, weights=T)
wAICc.body
summary(mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.body)

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

mod.fc_ed_mnnd.bill<-lmer(Bill~scale(fc)+scale(ed)+scale(mnnd)+ coords+ (1|ID_study),
                          data = dados, na.action = "na.fail")

####Calculating semivariogram

E <- resid(mod.fc_ed_mnnd.bill)
mydata <- data.frame(E, dados$lon, dados$lat)
coordinates(mydata) <- c("dados.lon", "dados.lat")

###semivariograma####
Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Bill width", cex=1.5, col="black")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.bill<-lmer(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                                data = dados, na.action = "na.fail")
summary(mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.bill)
mod.fc_ed_mnnd_int.fc_ed.bill<-lmer(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords+ (1|ID_study),
                                    data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.bill<-lmer(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                      data = dados, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.bill<-lmer(Bill~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study),
                                      data = dados, na.action = "na.fail")

mod.fc_ed.bill<-lmer(Bill~scale(fc)+scale(ed)+ coords+ (1|ID_study),
                     data = dados, na.action = "na.fail")

mod.fc_ed_int.fc_ed.bill<-lmer(Bill~scale(fc)+scale(ed)+scale(fc):scale(ed)+ coords+ (1|ID_study),
                               data = dados, na.action = "na.fail")

mod.fc_mnnd.bill<-lmer(Bill~scale(fc)+scale(mnnd)+ coords+ (1|ID_study),
                       data = dados, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.bill<-lmer(Bill~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                    data = dados, na.action = "na.fail")

mod.ed_mnnd.bill<-lmer(Bill~scale(ed)+scale(mnnd)+ coords+ (1|ID_study),
                       data = dados, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.bill<-lmer(Bill~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study),
                                   data = dados, na.action = "na.fail")

mod.fc.bill<-lmer(Bill~scale(fc)+ coords+ (1|ID_study),
                  data = dados, na.action = "na.fail")

mod.ed.bill<-lmer(Bill~scale(ed)+ coords+ (1|ID_study),
                  data = dados, na.action = "na.fail")

mod.mnnd.bill<-lmer(Bill~scale(mnnd)+ coords+ (1|ID_study),
                    data = dados, na.action = "na.fail")


wAICc.bill <- AICctab(mod.fc_ed_mnnd.bill, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.bill, mod.fc_ed_mnnd_int.fc_ed.bill,
                      mod.fc_ed_mnnd_int.fc_mnnd.bill, mod.fc_ed_mnnd_int.mnnd_ed.bill, mod.fc_ed.bill,
                      mod.fc_ed_int.fc_ed.bill, mod.fc_mnnd.bill, mod.fc_mnnd_int.fc_mnnnd.bill,
                      mod.ed_mnnd.bill,mod.ed_mnnd_int.ed_mnnd.bill, mod.fc.bill, mod.ed.bill, mod.mnnd.bill, weights=T)
wAICc.bill
summary(mod.ed_mnnd_int.ed_mnnd.bill)


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
str(dados2$Hand)
##For hand wing index we will exclude one line that has NA.
dados2<-dados[-c(130),]

dados2$coords<-as.matrix(cbind(dados2$lon,dados2$lat))

mod.fc_ed_mnnd.hand<-lmer(Hand~scale(fc)+scale(ed)+scale(mnnd)+ coords+ (1|ID_study),
                          data = dados2, na.action = "na.fail")

####Calculating semivariogram

E <- resid(mod.fc_ed_mnnd.hand)
mydata <- data.frame(E, dados2$lon, dados2$lat)
coordinates(mydata) <- c("dados2.lon", "dados2.lat")

###semivariograma####
Vario2 <- variogram(E ~ 1, mydata)
plot(Vario2, smooth = TRUE,pch=19, main = "Hand wing index", cex=1.5, col="black")

###Lets construct models
mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.hand<-lmer(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                                data = dados2, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_ed.hand<-lmer(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(ed)+ coords+ (1|ID_study),
                                    data = dados2, na.action = "na.fail")

mod.fc_ed_mnnd_int.fc_mnnd.hand<-lmer(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                      data = dados2, na.action = "na.fail")

mod.fc_ed_mnnd_int.mnnd_ed.hand<-lmer(Hand~scale(fc)+scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study),
                                      data = dados2, na.action = "na.fail")

mod.fc_ed.hand<-lmer(Hand~scale(fc)+scale(ed)+ coords+ (1|ID_study),
                     data = dados2, na.action = "na.fail")

mod.fc_ed_int.fc_ed.hand<-lmer(Hand~scale(fc)+scale(ed)+scale(fc):scale(ed)+ coords+ (1|ID_study),
                               data = dados2, na.action = "na.fail")

mod.fc_mnnd.hand<-lmer(Hand~scale(fc)+scale(mnnd)+ coords+ (1|ID_study),
                       data = dados2, na.action = "na.fail")

mod.fc_mnnd_int.fc_mnnnd.hand<-lmer(Hand~scale(fc)+scale(mnnd)+scale(fc):scale(mnnd)+ coords+ (1|ID_study),
                                    data = dados2, na.action = "na.fail")

mod.ed_mnnd.hand<-lmer(Hand~scale(ed)+scale(mnnd)+ coords+ (1|ID_study),
                       data = dados2, na.action = "na.fail")

mod.ed_mnnd_int.ed_mnnd.hand<-lmer(Hand~scale(ed)+scale(mnnd)+ scale(mnnd):scale(ed)+ coords+ (1|ID_study),
                                   data = dados2, na.action = "na.fail")

mod.fc.hand<-lmer(Hand~scale(fc)+ coords+ (1|ID_study),
                  data = dados2, na.action = "na.fail")

mod.ed.hand<-lmer(Hand~scale(ed)+ coords+ (1|ID_study),
                  data = dados2, na.action = "na.fail")

mod.mnnd.hand<-lmer(Hand~scale(mnnd)+ coords+ (1|ID_study),
                    data = dados2, na.action = "na.fail")


wAICc.hand <- AICctab(mod.fc_ed_mnnd.hand, mod.fc_ed_mnnd_int.fc_ed_int.fc_mnnd.hand, mod.fc_ed_mnnd_int.fc_ed.hand,
                      mod.fc_ed_mnnd_int.fc_mnnd.hand, mod.fc_ed_mnnd_int.mnnd_ed.hand, mod.fc_ed.hand,
                      mod.fc_ed_int.fc_ed.hand, mod.fc_mnnd.hand, mod.fc_mnnd_int.fc_mnnnd.hand,
                      mod.ed_mnnd.hand,mod.ed_mnnd_int.ed_mnnd.hand, mod.fc.hand, mod.ed.hand, mod.mnnd.hand, weights=T)
wAICc.hand

summary(mod.mnnd.hand)

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



