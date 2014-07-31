
rm(list=ls())

require(ggplot2)

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
source('collate.R')

# Load the Fishbase data
fb <- FishbaseWeb$read('data/fishbase-web')
# Load Gislason data
gs <- GislasonEtAl2010Data$read('data/gislason-et-al-2010')

# Create a fishnet that can be used to impute columns that 
# are missing in the Gislason data. This fishnet is intended to be 
# simple and mostly uses TaxonomicImputer
imputer <- Fishnet(
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit'),
  depthmax  = TaxonomicImputer('depthmax',c(log,exp)),
  temp      = TaxonomicImputer('temp'),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp))
)
# Do imputation base on Fishbase data
imputer$fit(fb)
gs <- imputer$predict(gs)

################
# Charnov 2013 #
################
source('m-charnov-et-al-2013.R')
source('m-charnov-et-al-2013-fitted.R')

# fishbase
cea13fb <- MCharnovEtAl2013()
fb.cea13fb <- data.frame(fb,m.predict=cea13fb$predict(fb))

cea13fbfit <- MCharnovEtAl2013Fitted()
cea13fbfit$fit(fb)
fb.cea13fbfit <- data.frame(fb,m.predict=cea13fbfit$predict(fb))

plot(log(m/k)~log(lmat/linf),data=fb,main='cea13fb')
points(log(m.predict/k)~log(lmat/linf),data=fb.cea13fb,col=2)
points(log(m.predict/k)~log(lmat/linf),data=fb.cea13fbfit,col=4)

cv.cea13fb <- cea13fb$cross(fb)
cv.cea13fbfit <- cea13fbfit$cross(fb)

# gislasson
gs$lmat = gs$l

cea13gs <- MCharnovEtAl2013()
gs.cea13gs <- data.frame(gs,m.predict=cea13gs$predict(gs))

cea13gsfit <- MCharnovEtAl2013Fitted()
cea13gsfit$fit(gs)
gs.cea13gsfit <- data.frame(gs,m.predict=cea13gsfit$predict(gs))

plot(log(m/k)~log(lmat/linf),data=gs,main='cea13gs')
points(log(m.predict/k)~log(lmat/linf),data=gs.cea13gs,col=2)
points(log(m.predict/k)~log(lmat/linf),data=gs.cea13gsfit,col=4)

cv.cea13gs <- cea13gs$cross(gs)
cv.cea13gsfit <- cea13gsfit$cross(gs)

cv.cea13 <- rbind(
  data.frame(source='cea13',method='GLM',fitted='N',db='fb',mpe=round(cv.cea13fb$summary$mpe   ,2),r2=round(cv.cea13fb$summary$r2,2)),
  data.frame(source='cea13',method='GLM',fitted='Y',db='fb',mpe=round(cv.cea13fbfit$summary$mpe,2),r2=round(cv.cea13fbfit$summary$r2,2)),
  data.frame(source='cea13',method='GLM',fitted='N',db='gs',mpe=round(cv.cea13gs$summary$mpe   ,2),r2=round(cv.cea13gs$summary$r2,2)),
  data.frame(source='cea13',method='GLM',fitted='Y',db='gs',mpe=round(cv.cea13gsfit$summary$mpe,2),r2=round(cv.cea13gsfit$summary$r2,2))
)


###############
# Hoenig 1983 #
###############
source('m-hoenig-1983.R')
source('m-hoenig-1983-fitted.R')

# fishbase
h83fb <- MHoenig1983()
fb.h83fb <- data.frame(fb,m.predict=h83fb$predict(fb))

h83fbfit <- MHoenig1983Fitted()
h83fbfit$fit(fb)
fb.h83fbfit <- data.frame(fb,m.predict=h83fbfit$predict(fb))

plot(log(m)~log(amax),data=fb,main='h83fb')
points(log(m.predict)~log(amax),data=fb.h83fb,col=2)
points(log(m.predict)~log(amax),data=fb.h83fbfit,col=4)

cv.h83fb <- h83fb$cross(fb)
cv.h83fbfit <- h83fbfit$cross(fb)

# gislasson

h83gs <- MHoenig1983()
gs.h83gs <- data.frame(gs,m.predict=h83gs$predict(gs))

h83gsfit <- MHoenig1983Fitted()
h83gsfit$fit(gs)
gs.h83gsfit <- data.frame(gs,m.predict=h83gsfit$predict(gs))

plot(log(m)~log(amax),data=gs,main='h83gs')
points(log(m.predict)~log(amax),data=gs.h83gs,col=2)
points(log(m.predict)~log(amax),data=gs.h83gsfit,col=4)

cv.h83gs <- h83gs$cross(gs)
cv.h83gsfit <- h83gsfit$cross(gs)

cv.h83 <- rbind(
  data.frame(source='h83',method='GLM',fitted='N',db='fb',mpe=round(cv.h83fb$summary$mpe   ,2),r2=round(cv.h83fb$summary$r2,2)),
  data.frame(source='h83',method='GLM',fitted='Y',db='fb',mpe=round(cv.h83fbfit$summary$mpe,2),r2=round(cv.h83fbfit$summary$r2,2)),
  data.frame(source='h83',method='GLM',fitted='N',db='gs',mpe=round(cv.h83gs$summary$mpe   ,2),r2=round(cv.h83gs$summary$r2,2)),
  data.frame(source='h83',method='GLM',fitted='Y',db='gs',mpe=round(cv.h83gsfit$summary$mpe,2),r2=round(cv.h83gsfit$summary$r2,2))
)

#############
# Then 2014 #
#############

source('m-then-et-al-2014.R')
source('m-then-et-al-2014-fitted.R')

# fishbase
tea14fb <- MThenEtAl2014()
fb.tea14fb <- data.frame(fb,m.predict=tea14fb$predict(fb))

tea14fbfit <- MThenEtAl2014Fitted()
tea14fbfit$fit(fb)
fb.tea14fbfit <- data.frame(fb,m.predict=tea14fbfit$predict(fb))

plot(log(m)~log(amax),data=fb,main='tea14fb')
points(log(m.predict)~log(amax),data=fb.tea14fbfit,col=4)

cv.tea14fbfit <- tea14fbfit$cross(fb)

# gislasson

tea14gs <- MHoenig1983()
gs.tea14gs <- data.frame(gs,m.predict=tea14gs$predict(gs))

tea14gsfit <- MHoenig1983Fitted()
tea14gsfit$fit(gs)
gs.tea14gsfit <- data.frame(gs,m.predict=tea14gsfit$predict(gs))

plot(log(m)~log(amax),data=gs,main='tea14gs')
points(log(m.predict)~log(amax),data=gs.tea14gsfit,col=4)

cv.tea14gsfit <- tea14gsfit$cross(gs)

cv.tea14 <- rbind(
  data.frame(source='tea14',method='GLM',fitted='N',db='fb',mpe=NA                                ,r2=NA),
  data.frame(source='tea14',method='GLM',fitted='Y',db='fb',mpe=round(cv.tea14fbfit$summary$mpe,2),r2=round(cv.tea14fbfit$summary$r2,2)),
  data.frame(source='tea14',method='GLM',fitted='N',db='gs',mpe=NA                                ,r2=NA),
  data.frame(source='tea14',method='GLM',fitted='Y',db='gs',mpe=round(cv.tea14gsfit$summary$mpe,2),r2=round(cv.tea14gsfit$summary$r2,2))
)

##############
# Pauly 1980 #
##############
source('m-pauly-1980.R')
source('m-pauly-1980-fitted.R')

# fishbase
p80fb <- MPauly1980()
fb.p80fb <- data.frame(fb,m.predict=p80fb$predict(fb))

p80fbfit <- MPauly1980Fitted()
p80fbfit$fit(fb)
fb.p80fbfit <- data.frame(fb,m.predict=p80fbfit$predict(fb))

plot(log(m/k)~log(linf),data=fb,main='p80fb')
points(log(m.predict/k)~log(linf),data=fb.p80fb,col=2)
points(log(m.predict/k)~log(linf),data=fb.p80fbfit,col=4)

cv.p80fb <- p80fb$cross(fb)
cv.p80fbfit <- p80fbfit$cross(fb)

# gislasson

p80gs <- MPauly1980()
gs.p80gs <- data.frame(gs,m.predict=p80gs$predict(gs))

p80gsfit <- MPauly1980Fitted()
p80gsfit$fit(gs)
gs.p80gsfit <- data.frame(gs,m.predict=p80gsfit$predict(gs))

plot(log(m/k)~log(linf),data=gs,main='p80gs')
points(log(m.predict/k)~log(linf),data=gs.p80gs,col=2)
points(log(m.predict/k)~log(linf),data=gs.p80gsfit,col=4)

cv.p80gs <- p80gs$cross(gs)
cv.p80gsfit <- p80gsfit$cross(gs)

cv.p80 <- rbind(
  data.frame(source='p80',method='GLM',fitted='N',db='fb',mpe=round(cv.p80fb$summary$mpe   ,2),r2=round(cv.p80fb$summary$r2,2)),
  data.frame(source='p80',method='GLM',fitted='Y',db='fb',mpe=round(cv.p80fbfit$summary$mpe,2),r2=round(cv.p80fbfit$summary$r2,2)),
  data.frame(source='p80',method='GLM',fitted='N',db='gs',mpe=round(cv.p80gs$summary$mpe   ,2),r2=round(cv.p80gs$summary$r2,2)),
  data.frame(source='p80',method='GLM',fitted='Y',db='gs',mpe=round(cv.p80gsfit$summary$mpe,2),r2=round(cv.p80gsfit$summary$r2,2))
)



###########
# SUMMARY #
###########

cv.sum <- rbind(cv.cea13,cv.h83,cv.tea14,cv.p80)


# boosted regression trees
brt.min <- Brter(log(m)~log(k)+log(amax),exp)
fb.min <- subset(fb,!is.na(k) & !is.na(amax))
brt.min$fit(fb.min)
brt.min$cross(fb.min)
hist(brt.min$sample(brt.min$expand(data.frame(k=.3,amax=5),1000)))
summary(brt.min$brt)
nrow(fb.min)

brt.full <- Brter(log(m)~class+order+family+log(linf)+log(temp)+log(k)+log(amax),exp)
fb.full <- subset(fb,!is.na(class) & !is.na(order) & !is.na(family) & !is.na(linf) & !is.na(temp) & temp>0 & !is.na(k) & !is.na(amax))
brt.full$fit(fb.full)
brt.full$cross(fb.full)
hist(brt.full$sample(brt.full$expand(data.frame(k=.3,amax=5),1000)))
summary(brt.full$brt)
nrow(fb.full)
