
rm(list=ls())

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
source('collate.R')
source('load_data.R')

################
# Charnov 2013 #
################

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
  data.frame(source='cea13',method='GLM',fitted=FALSE,db='fb',mpe=round(cv.cea13fb$summary['mpe',1]   ,2),dev=round(cv.cea13fb$summary['dev',1],2)),
  data.frame(source='cea13',method='GLM',fitted=TRUE,db='fb',mpe=round(cv.cea13fbfit$summary['mpe',1],2),dev=round(cv.cea13fbfit$summary['dev',1],2)),
  data.frame(source='cea13',method='GLM',fitted=FALSE,db='gs',mpe=round(cv.cea13gs$summary['mpe',1]   ,2),dev=round(cv.cea13gs$summary['dev',1],2)),
  data.frame(source='cea13',method='GLM',fitted=TRUE,db='gs',mpe=round(cv.cea13gsfit$summary['mpe',1],2),dev=round(cv.cea13gsfit$summary['dev',1],2))
)


###############
# Hoenig 1983 #
###############

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
  data.frame(source='h83',method='GLM',fitted=FALSE,db='fb',mpe=round(cv.h83fb$summary['mpe',1]   ,2),dev=round(cv.h83fb$summary['dev',1],2)),
  data.frame(source='h83',method='GLM',fitted=TRUE,db='fb',mpe=round(cv.h83fbfit$summary['mpe',1],2),dev=round(cv.h83fbfit$summary['dev',1],2)),
  data.frame(source='h83',method='GLM',fitted=FALSE,db='gs',mpe=round(cv.h83gs$summary['mpe',1]   ,2),dev=round(cv.h83gs$summary['dev',1],2)),
  data.frame(source='h83',method='GLM',fitted=TRUE,db='gs',mpe=round(cv.h83gsfit$summary['mpe',1],2),dev=round(cv.h83gsfit$summary['dev',1],2))
)

#############
# Then 2014 #
#############

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

tea14gs <- MThenEtAl2014()
gs.tea14gs <- data.frame(gs,m.predict=tea14gs$predict(gs))

tea14gsfit <- MThenEtAl2014Fitted()
tea14gsfit$fit(gs)
gs.tea14gsfit <- data.frame(gs,m.predict=tea14gsfit$predict(gs))

plot(log(m)~log(amax),data=gs,main='tea14gs')
points(log(m.predict)~log(amax),data=gs.tea14gsfit,col=4)

cv.tea14gsfit <- tea14gsfit$cross(gs)

cv.tea14 <- rbind(
  data.frame(source='tea14',method='GLM',fitted=FALSE,db='fb',mpe=NA                                ,dev=NA),
  data.frame(source='tea14',method='GLM',fitted=TRUE,db='fb',mpe=round(cv.tea14fbfit$summary['mpe',1],2),dev=round(cv.tea14fbfit$summary['dev',1],2)),
  data.frame(source='tea14',method='GLM',fitted=FALSE,db='gs',mpe=NA                                ,dev=NA),
  data.frame(source='tea14',method='GLM',fitted=TRUE,db='gs',mpe=round(cv.tea14gsfit$summary['mpe',1],2),dev=round(cv.tea14gsfit$summary['dev',1],2))
)

##############
# Pauly 1980 #
##############

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
  data.frame(source='p80',method='GLM',fitted=FALSE,db='fb',mpe=round(cv.p80fb$summary['mpe',1]   ,2),dev=round(cv.p80fb$summary['dev',1],2)),
  data.frame(source='p80',method='GLM',fitted=TRUE,db='fb',mpe=round(cv.p80fbfit$summary['mpe',1],2),dev=round(cv.p80fbfit$summary['dev',1],2)),
  data.frame(source='p80',method='GLM',fitted=FALSE,db='gs',mpe=round(cv.p80gs$summary['mpe',1]   ,2),dev=round(cv.p80gs$summary['dev',1],2)),
  data.frame(source='p80',method='GLM',fitted=TRUE,db='gs',mpe=round(cv.p80gsfit$summary['mpe',1],2),dev=round(cv.p80gsfit$summary['dev',1],2))
)

###########
# SUMMARY #
###########

cv.sum <- rbind(cv.cea13,cv.h83,cv.tea14,cv.p80)
subset(cv.sum,fitted=='Y')


############################
# boosted regression trees #
############################

# fb

# model selection using gbm.simplify()
frame <- model.frame(log(m)~class+order+family+log(linf)+log(lmat)+log(temp)+log(k)+log(amax),fb)
names(frame) <- c("log.m","class","order","family","log.linf","log.lmat","log.temp","log.k","log.amax")
brtfb <- gbm.step(data = frame,
                  gbm.y = 1,
                  gbm.x = 2:ncol(frame), 
                  family = "gaussian",
                  tree.complexity = 10,
                  learning.rate = 0.002,
                  bag.fraction = 0.5,
                  max.trees = 5000)

summary(brtfb)

brtfb.simplify <- gbm.simplify(brtfb,n.drops=length(brtfb$gbm.call$gbm.x)-2)

brtfb.simplify$final.drops

# let's try within fishnets

brtfb <- Brter(formula=log(m)~class+order+family+log(linf)+log(lmat)+log(temp)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)

brtfb$fit(fb)
pars <- brtfb$brt$gbm.call$gbm.x

res <- list()
pars.update <- pars

for(i in 1:(length(pars)-1)) {
  
 brtfb$fit(fb,pars.update)
 cat('predictors:',brtfb$brt$gbm.call$predictor.names,'\n')
 
 rinfl <- relative.influence(brtfb$brt,n.trees = brtfb$brt$gbm.call$best.trees)
 
 res[[i]] <- list()
 res[[i]][['summary']]   <- brtfb$cross(fb,pars=pars.update)$summary 
 res[[i]][['influence']] <- data.frame(predictor=brtfb$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
 res[[i]][['drop']]      <- brtfb$brt$gbm.call$predictor.names[which.min(rinfl)]
 
 
 pars.update <- pars.update[-which.min(rinfl)]
 
}

# results are consistent with gbm.simplify()
unlist(lapply(res,function(x) x$drop)

# ok now manually... we need to do this to use all the data with each drop (sigh!)

# full
brtfb <- Brter(formula=log(m)~class+order+family+log(linf)+log(lmat)+log(temp)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[1]][['summary2']]   <- brtfb$cross(fb)$summary 
# drop class
brtfb <- Brter(formula=log(m)~order+family+log(linf)+log(lmat)+log(temp)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[2]][['summary2']]   <- brtfb$cross(fb)$summary 
# drop log(temp)
brtfb <- Brter(formula=log(m)~order+family+log(linf)+log(lmat)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[3]][['summary2']]   <- brtfb$cross(fb)$summary 
# drop order
brtfb <- Brter(formula=log(m)~family+log(linf)+log(lmat)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[4]][['summary2']]   <- brtfb$cross(fb)$summary 
# drop log(lmat)
brtfb <- Brter(formula=log(m)~family+log(linf)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[5]][['summary2']]   <- brtfb$cross(fb)$summary 
# drop log(linf)
brtfb <- Brter(formula=log(m)~family+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[6]][['summary2']]   <- brtfb$cross(fb)$summary 
# drop log(k)
brtfb <- Brter(formula=log(m)~family+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
res[[7]][['summary2']]   <- brtfb$cross(fb)$summary 

dev <- unlist(lapply(res,function(x) x$summary2['dev','mean']))
dev.se <- unlist(lapply(res,function(x) x$summary2['dev','se']))

plot(1:length(predictor),dev,type='l',xaxt='n',ylim=range(dev-1.1*dev.se,dev+1.1*dev.se),xlab='Predictor removed',ylab='Deviance')
lines(1:length(predictor),dev+dev.se,lty=2)
lines(1:length(predictor),dev-dev.se,lty=2)
abline(h=res[[1]]$summary2['dev','mean'],lty=2,col=2)
axis(1,at=1:length(predictor),labels=c('none',predictor[-length(predictor)]))

# FINAL MODEL
brtfb <- Brter(formula=log(m)~family+log(linf)+log(lmat)+log(k)+log(amax),transform = exp,ntrees=0,learning.rate=0.002,max.trees=10000)
brtfb$fit(fb)
save(brtfb.final=brtfb,res,file='examples/m-trials-cvresults.Rdata')










