
#######################################
# DETERMINISTIC RELATIONSHIPS BETWEEN #
# M AND OTHER LH PARMETERS            #
#######################################

rm(list=ls())

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
source('collate.R')
source('load_data.R')
source('utils.R')

##############
# check data #
##############

# counts
fb.tmp <- fb[!is.na(fb$m),]
for(i in 1:ncol(fb)) {
  print(paste(names(fb)[i],':',sum(!is.na(fb.tmp[,i]))))
}
rm(fb.tmp)

gs.tmp <- gs[!is.na(gs$m),]
for(i in 1:ncol(gs)) {
  print(paste(names(gs)[i],':',sum(!is.na(gs.tmp[,i]))))
}
rm(gs.tmp)

# check that taxanomic imputation of lmat corresponds
# to that expected using BH invariant in gislasson database
# (Charnov 2013)
plot(gs$lmat,gs$linf * 2/3); abline(0,1)
gs$lmat[gs$lmat>150] <- NA

################
# Charnov 1990 #
################

# fishbase
c90fb <- MCharnov1990()
fb.c90fb <- data.frame(fb,m.predict=c90fb$predict(fb))

c90fbfit <- MCharnov1990Fitted()
c90fbfit$fit(fb)
fb.c90fbfit <- data.frame(fb,m.predict=c90fbfit$predict(fb))

#plot(log(m)~log(amat),data=fb,main='c90fb')
#points(log(m.predict)~log(amat),data=fb.c90fb,col=2)
#points(log(m.predict)~log(amat),data=fb.c90fbfit,col=4)

cv.c90fb    <- c90fb$cross(fb)
cv.c90fbfit <- c90fbfit$cross(fb)

# gislason
c90gs <- MCharnov1990()
gs.c90gs <- data.frame(gs,m.predict=c90gs$predict(gs))

c90gsfit <- MCharnov1990Fitted()
c90gsfit$fit(gs)
gs.c90gsfit <- data.frame(gs,m.predict=c90gsfit$predict(gs))

par(mfrow=c(1,2),mar=c(3,2,1,1))
plot(log(m)~log(amat),data=gs,main='c90gs')
points(log(m.predict)~log(amat),data=gs.c90gs,col=2)
points(log(m.predict)~log(amat),data=gs.c90gsfit,col=4)
plot(m.predict~m,data=gs.c90gsfit,main='c90gs')
abline(0,1,col=2)

cv.c90gs    <- c90gs$cross(gs)
cv.c90gsfit <- c90gsfit$cross(gs)

# summary
cv.c90 <- rbind(
  data.frame(source='(1)',method='GLM',fitted=FALSE,db='fb',n=NA,mpe=round(cv.c90fb$summary['mpe',1],2),dev=round(cv.c90fb$summary['dev',1],2)),
  data.frame(source='(1)',method='GLM',fitted=TRUE,db='fb',n=c90fbfit$n(fb),mpe=round(cv.c90fbfit$summary['mpe',1],2),dev=round(cv.c90fbfit$summary['dev',1],2)),
  data.frame(source='(1)',method='GLM',fitted=FALSE,db='gs',n=NA,mpe=round(cv.c90gs$summary['mpe',1],2),dev=round(cv.c90gs$summary['dev',1],2)),
  data.frame(source='(1)',method='GLM',fitted=TRUE,db='gs',n=c90gsfit$n(gs),mpe=round(cv.c90gsfit$summary['mpe',1],2),dev=round(cv.c90gsfit$summary['dev',1],2))
)

################
# Charnov 1993 #
################

# fishbase
c93fb <- MCharnov1993()
fb.c93fb <- data.frame(fb,m.predict=c93fb$predict(fb))

c93fbfit <- MCharnov1993Fitted()
c93fbfit$fit(fb)
fb.c93fbfit <- data.frame(fb,m.predict=c93fbfit$predict(fb))

#plot(log(m)~log(k),data=fb,main='c93fb')
#points(log(m.predict)~log(k),data=fb.c93fb,col=2)
#points(log(m.predict)~log(k),data=fb.c93fbfit,col=4)

cv.c93fb    <- c93fb$cross(fb)
cv.c93fbfit <- c93fbfit$cross(fb)

# gislason
c93gs <- MCharnov1993()
gs.c93gs <- data.frame(gs,m.predict=c93gs$predict(gs))

c93gsfit <- MCharnov1993Fitted()
c93gsfit$fit(gs)
gs.c93gsfit <- data.frame(gs,m.predict=c93gsfit$predict(gs))

#plot(log(m)~log(k),data=gs,main='c93gs')
#points(log(m.predict)~log(k),data=gs.c93gs,col=2)
#points(log(m.predict)~log(k),data=gs.c93gsfit,col=4)

par(mfrow=c(1,2),mar=c(3,2,1,1))
plot(log(m)~log(k),data=gs,main='c93gs')
points(log(m.predict)~log(k),data=gs.c93gs,col=2)
points(log(m.predict)~log(k),data=gs.c93gsfit,col=4)
plot(m.predict~m,data=gs.c93gsfit,main='c93gs')
abline(0,1,col=2)

cv.c93gs    <- c93gs$cross(gs)
cv.c93gsfit <- c93gsfit$cross(gs)

# summary
cv.c93 <- rbind(
  data.frame(source='(2)',method='GLM',fitted=FALSE,db='fb',n=NA,mpe=round(cv.c93fb$summary['mpe',1],2),dev=round(cv.c93fb$summary['dev',1],2)),
  data.frame(source='(2)',method='GLM',fitted=TRUE,db='fb',n=c93fbfit$n(fb),mpe=round(cv.c93fbfit$summary['mpe',1],2),dev=round(cv.c93fbfit$summary['dev',1],2)),
  data.frame(source='(2)',method='GLM',fitted=FALSE,db='gs',n=NA,mpe=round(cv.c93gs$summary['mpe',1],2),dev=round(cv.c93gs$summary['dev',1],2)),
  data.frame(source='(2)',method='GLM',fitted=TRUE,db='gs',n=c93gsfit$n(gs),mpe=round(cv.c93gsfit$summary['mpe',1],2),dev=round(cv.c93gsfit$summary['dev',1],2))
)


################
# Charnov 2013 #
################

# fishbase
c13fb <- MCharnov2013()
fb.c13fb <- data.frame(fb,m.predict=c13fb$predict(fb))

c13fbfit <- MCharnov2013Fitted()
c13fbfit$fit(fb)
fb.c13fbfit <- data.frame(fb,m.predict=c13fbfit$predict(fb))

#plot(log(m/k)~log(lmat/linf),data=fb,main='c13fb')
#points(log(m.predict/k)~log(lmat/linf),data=fb.c13fb,col=2)
#points(log(m.predict/k)~log(lmat/linf),data=fb.c13fbfit,col=4)

cv.c13fb <- c13fb$cross(fb)
cv.c13fbfit <- c13fbfit$cross(fb)

# gislason
c13gs <- MCharnov2013()
gs.c13gs <- data.frame(gs,m.predict=c13gs$predict(gs))

c13gsfit <- MCharnov2013Fitted()
c13gsfit$fit(gs)
gs.c13gsfit <- data.frame(gs,m.predict=c13gsfit$predict(gs))

par(mfrow=c(1,2),mar=c(3,2,1,1))
plot(log(m/k)~log(lmat/linf),data=gs,main='c13gs')
points(log(m.predict/k)~log(lmat/linf),data=gs.c13gs,col=2)
points(log(m.predict/k)~log(lmat/linf),data=gs.c13gsfit,col=4)
plot(m.predict~m,data=gs.c13gsfit,main='c13gs')
abline(0,1,col=2)


cv.c13gs <- c13gs$cross(gs)
cv.c13gsfit <- c13gsfit$cross(gs)

# summary
cv.c13 <- rbind(
  data.frame(source='(5)',method='GLM',fitted=FALSE,db='fb',n=NA,mpe=round(cv.c13fb$summary['mpe',1],2),dev=round(cv.c13fb$summary['dev',1],2)),
  data.frame(source='(5)',method='GLM',fitted=TRUE,db='fb',n=c13fbfit$n(fb),mpe=round(cv.c13fbfit$summary['mpe',1],2),dev=round(cv.c13fbfit$summary['dev',1],2)),
  data.frame(source='(5)',method='GLM',fitted=FALSE,db='gs',n=NA,mpe=round(cv.c13gs$summary['mpe',1],2),dev=round(cv.c13gs$summary['dev',1],2)),
  data.frame(source='(5)',method='GLM',fitted=TRUE,db='gs',n=c13gsfit$n(gs),mpe=round(cv.c13gsfit$summary['mpe',1],2),dev=round(cv.c13gsfit$summary['dev',1],2))
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

#plot(log(m)~log(amax),data=fb,main='h83fb')
#points(log(m.predict)~log(amax),data=fb.h83fb,col=2)
#points(log(m.predict)~log(amax),data=fb.h83fbfit,col=4)

cv.h83fb <- h83fb$cross(fb)
cv.h83fbfit <- h83fbfit$cross(fb)

# gislasson

h83gs <- MHoenig1983()
gs.h83gs <- data.frame(gs,m.predict=h83gs$predict(gs))

h83gsfit <- MHoenig1983Fitted()
h83gsfit$fit(gs)
gs.h83gsfit <- data.frame(gs,m.predict=h83gsfit$predict(gs))

#plot(log(m)~log(amax),data=gs,main='h83gs')
#points(log(m.predict)~log(amax),data=gs.h83gs,col=2)
#points(log(m.predict)~log(amax),data=gs.h83gsfit,col=4)

cv.h83gs <- h83gs$cross(gs)
cv.h83gsfit <- h83gsfit$cross(gs)

cv.h83 <- rbind(
  data.frame(source='(3)',method='GLM',fitted=FALSE,db='fb',n=NA,mpe=round(cv.h83fb$summary['mpe',1]   ,2),dev=round(cv.h83fb$summary['dev',1],2)),
  data.frame(source='(3)',method='GLM',fitted=TRUE,db='fb',n=h83fbfit$n(fb),mpe=round(cv.h83fbfit$summary['mpe',1],2),dev=round(cv.h83fbfit$summary['dev',1],2)),
  data.frame(source='(3)',method='GLM',fitted=FALSE,db='gs',n=NA,mpe=round(cv.h83gs$summary['mpe',1]   ,2),dev=round(cv.h83gs$summary['dev',1],2)),
  data.frame(source='(3)',method='GLM',fitted=TRUE,db='gs',n=h83gsfit$n(gs),mpe=round(cv.h83gsfit$summary['mpe',1],2),dev=round(cv.h83gsfit$summary['dev',1],2))
)

##############
# Quinn 1999 #
##############

# fishbase
q99fb <- MQuinn1999()
fb.q99fb <- data.frame(fb,m.predict=q99fb$predict(fb))

q99fbfit <- MQuinn1999Fitted()
q99fbfit$fit(fb)
fb.q99fbfit <- data.frame(fb,m.predict=q99fbfit$predict(fb))

#plot(log(m)~log(amax),data=fb,main='q99fb')
#points(log(m.predict)~log(amax),data=fb.q99fb,col=2)
#points(log(m.predict)~log(amax),data=fb.q99fbfit,col=4)

cv.q99fb <- q99fb$cross(fb)
cv.q99fbfit <- q99fbfit$cross(fb)

# gislasson
q99gs <- MQuinn1999()
gs.q99gs <- data.frame(gs,m.predict=q99gs$predict(gs))

q99gsfit <- MQuinn1999Fitted()
q99gsfit$fit(gs)
gs.q99gsfit <- data.frame(gs,m.predict=q99gsfit$predict(gs))

#plot(log(m)~log(amax),data=gs,main='q99gs')
#points(log(m.predict)~log(amax),data=gs.q99gs,col=2)
#points(log(m.predict)~log(amax),data=gs.q99gsfit,col=4)

cv.q99gs <- q99gs$cross(gs)
cv.q99gsfit <- q99gsfit$cross(gs)

cv.q99 <- rbind(
  data.frame(source='(4)',method='GLM',fitted=FALSE,db='fb',n=NA,mpe=round(cv.q99fb$summary['mpe',1],2),dev=round(cv.q99fb$summary['dev',1],2)),
  data.frame(source='(4)',method='GLM',fitted=TRUE,db='fb',n=q99fbfit$n(fb),mpe=round(cv.q99fbfit$summary['mpe',1],2),dev=round(cv.q99fbfit$summary['dev',1],2)),
  data.frame(source='(4)',method='GLM',fitted=FALSE,db='gs',n=NA,mpe=round(cv.q99gs$summary['mpe',1],2),dev=round(cv.q99fb$summary['dev',1],2)),
  data.frame(source='(4)',method='GLM',fitted=TRUE,db='gs',n=q99gsfit$n(gs),mpe=round(cv.q99gsfit$summary['mpe',1],2),dev=round(cv.q99gsfit$summary['dev',1],2))
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

#plot(log(m/k)~log(linf),data=fb,main='p80fb')
#points(log(m.predict/k)~log(linf),data=fb.p80fb,col=2)
#points(log(m.predict/k)~log(linf),data=fb.p80fbfit,col=4)

cv.p80fb <- p80fb$cross(fb)
cv.p80fbfit <- p80fbfit$cross(fb)

# gislasson

p80gs <- MPauly1980()
gs.p80gs <- data.frame(gs,m.predict=p80gs$predict(gs))

p80gsfit <- MPauly1980Fitted()
p80gsfit$fit(gs)
gs.p80gsfit <- data.frame(gs,m.predict=p80gsfit$predict(gs))

cv.p80gs <- p80gs$cross(gs)
cv.p80gsfit <- p80gsfit$cross(gs)

cv.p80 <- rbind(
  data.frame(source='(6)',method='GLM',fitted=FALSE,db='fb',n=NA,mpe=round(cv.p80fb$summary['mpe',1],2),dev=round(cv.p80fb$summary['dev',1],2)),
  data.frame(source='(6)',method='GLM',fitted=TRUE,db='fb',n=p80fbfit$n(fb),mpe=round(cv.p80fbfit$summary['mpe',1],2),dev=round(cv.p80fbfit$summary['dev',1],2)),
  data.frame(source='(6)',method='GLM',fitted=FALSE,db='gs',n=NA,mpe=round(cv.p80gs$summary['mpe',1],2),dev=round(cv.p80gs$summary['dev',1],2)),
  data.frame(source='(6)',method='GLM',fitted=TRUE,db='gs',n=p80gsfit$n(gs),mpe=round(cv.p80gsfit$summary['mpe',1],2),dev=round(cv.p80gsfit$summary['dev',1],2))
)

###############
# Jensen 2001 #
###############

# gislasson

j01gs <- MJensen2001()
gs.j01gs <- data.frame(gs,m.predict=j01gs$predict(gs))

j01gsfit <- MJensen2001Fitted()
j01gsfit$fit(gs)
gs.j01gsfit <- data.frame(gs,m.predict=j01gsfit$predict(gs))

par(mfrow=c(1,2),mar=c(3,2,1,1))
plot(log(m)~log(k),data=gs,main='j01gs')
points(log(m.predict)~log(k),data=gs.j01gs,col=2)
points(log(m.predict)~log(k),data=gs.j01gsfit,col=4)
plot(m.predict~m,data=gs.j01gsfit,main='j01gs')
abline(0,1,col=2)


cv.j01gs <- j01gs$cross(gs)
cv.j01gsfit <- j01gsfit$cross(gs)


###########
# SUMMARY #
###########

# collate results
c90 <- list()
c90$fb <- list()
c90$fb$lit <- c90fb 
c90$fb$fit <- c90fbfit 
c90$gs <- list()
c90$gs$lit <- c90gs 
c90$gs$fit <- c90gsfit 

c93 <- list()
c93$fb <- list()
c93$fb$lit <- c93fb 
c93$fb$fit <- c93fbfit 
c93$gs <- list()
c93$gs$lit <- c93gs 
c93$gs$fit <- c93gsfit 

h83 <- list()
h83$fb <- list()
h83$fb$lit <- h83fb 
h83$fb$fit <- h83fbfit 
h83$gs <- list()
h83$gs$lit <- h83gs 
h83$gs$fit <- h83gsfit 

q99 <- list()
q99$fb <- list()
q99$fb$lit <- q99fb 
q99$fb$fit <- q99fbfit 
q99$gs <- list()
q99$gs$lit <- q99gs 
q99$gs$fit <- q99gsfit 

c13 <- list()
c13$fb <- list()
c13$fb$lit <- c13fb 
c13$fb$fit <- c13fbfit 
c13$gs <- list()
c13$gs$lit <- c13gs 
c13$gs$fit <- c13gsfit 

p80 <- list()
p80$fb <- list()
p80$fb$lit <- p80fb 
p80$fb$fit <- p80fbfit 
p80$gs <- list()
p80$gs$lit <- p80gs 
p80$gs$fit <- p80gsfit 

# cross validation tables
dfr <- rbind(cv.c90,cv.c93,cv.h83,cv.q99,cv.c13,cv.p80)
dfr <- dfr[,-2]
dfr <- rbind(subset(dfr,!fitted),subset(dfr,fitted))

# regression coefficients
exp(coef(c90fbfit$glm))
exp(coef(c93fbfit$glm))
exp(coef(h83fbfit$glm)[1]);coef(h83fbfit$glm)[2]
exp(coef(q99fbfit$glm))
coef(c13fbfit$glm)
coef(p80fbfit$glm)

exp(coef(c90gsfit$glm))
exp(coef(c93gsfit$glm))
exp(coef(h83gsfit$glm)[1]);coef(h83gsfit$glm)[2]
exp(coef(q99gsfit$glm))
coef(c13gsfit$glm)
coef(p80gsfit$glm)

# save
saver(c90,c93,h83,q99,c13,p80,dfr,name='m_det_res')

write.csv(dfr,file='C:/PROJECTS/FISHNETS/res/cvdet.csv')














