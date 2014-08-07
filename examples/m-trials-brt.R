
############################
# boosted regression trees #
############################

library(ggplot2)

rm(list=ls())
source('collate.R')
source('load_data.R')

saver <- function(x, ..., name, path='C:/PROJECTS/FISHNETS/res/') {
  save(x, ..., file=paste(path,name,'.Rdata',sep=''))
}

# according to life history theory, primary correlates
# with m are: k and amax; via the BH invariants
# m/k and m*amax. Therefore need to optimise estimation
# of these too.

##############
# groom data #
##############

fb[which(fb$temp<=0),'temp'] <- NA
fb[which(fb$m>2),'m']        <- NA

############################
# amax                     #
############################

formula <- log(amax)~class+order+family+sex+trophic+habit+log(temp)+log(lmax)
vars    <- all.vars(formula)
frame   <- model.frame(formula,fb)

names(frame) <- c("log.amax","class","order","family","sex","trophic","habit","log.temp","log.lmax")
brt.amax <- gbm.step(data = frame,
                  gbm.y = 1,
                  gbm.x = 2:ncol(frame), 
                  family = "gaussian",
                  tree.complexity = 10,
                  learning.rate = 0.004,
                  bag.fraction = 0.5,
                  max.trees = 10000)

summary(brt.amax)

brt.amax.simplify <- gbm.simplify(brt.amax)
lapply(brt.amax.simplify$pred.list,function(x) vars[x])

brt.amax.predictors <- vars[brt.amax.simplify$pred.list[[3]]]
brt.amax.formula    <- log(amax)~family+trophic+habit+log(temp)+log(lmax)

brt.amax <- Brter(brt.amax.formula,exp,tree.complexity = 10,learning.rate = 0.004,bag.fraction = 0.5,ntrees = 5900)
brt.amax$fit(fb)

brt.amax.cv <- brt.amax$cross(fb,folds=100)

dfr <- cbind(par='amax',brt.amax.cv$folds[,c('hat','obs')])

ggplot(dfr) + 
  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  theme_bw(base_size=20) + 
  labs(x='Observed',y='Predicted')

saver(brt.amax,brt.amax.formula,brt.amax.predictors,brt.amax.cv,name='brt_amax')

############################
# k                        #
############################

# NB this node performs slightly better if we use linf instead of lmax
# but looking at lmax and linf they are very closely related and better
# to use a variable that is measured directly rather than estimated

formula <- log(k)~class+order+family+sex+trophic+habit+log(temp)+log(lmax)+log(amax)
vars    <- all.vars(formula)
frame   <- model.frame(formula,fb)

names(frame) <- c("log.k","class","order","family","sex","trophic","habit","log.temp","log.lmax","log.amax")
brt.k <- gbm.step(data = frame,
                     gbm.y = 1,
                     gbm.x = 2:ncol(frame), 
                     family = "gaussian",
                     tree.complexity = 10,
                     learning.rate = 0.004,
                     bag.fraction = 0.5,
                     max.trees = 10000)

summary(brt.k)

brt.k.simplify <- gbm.simplify(brt.k)
lapply(brt.k.simplify$pred.list,function(x) vars[x])

brt.k.predictors <- vars[brt.k.simplify$pred.list[[4]]]
brt.k.formula    <- log(k)~family+trophic+log(temp)+log(lmax)+log(amax)

brt.k <- Brter(brt.k.formula,exp,tree.complexity = 10,learning.rate = 0.004,bag.fraction = 0.5,ntrees = 2400)
#brt.k$fit(fb)

brt.k.cv <- brt.k$cross(fb,folds=10)

dfr <- cbind(par='k',brt.k.cv$folds[,c('hat','obs')])

ggplot(dfr) + 
  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  theme_bw(base_size=20) + 
  labs(x='Observed',y='Predicted')

saver(brt.k,brt.k.formula,brt.k.predictors,brt.k.cv,name='brt_k')

############################
# m                        #
############################

# NB again this might work better with linf instead of lmax

formula <- log(m)~class+order+family+sex+trophic+habit+log(temp)+log(lmax)+log(k)+log(amax)
vars    <- all.vars(formula)
frame <- model.frame(formula,fb)
names(frame) <- c("log.m","class","order","family","sex","trophic","habit","log.temp","log.lmax","log.k","log.amax")
brt.m <- gbm.step(data = frame,
                  gbm.y = 1,
                  gbm.x = 2:ncol(frame), 
                  family = "gaussian",
                  tree.complexity = 10,
                  learning.rate = 0.002,
                  bag.fraction = 0.5,
                  max.trees = 10000)

summary(brt.m)

brt.m.simplify <- gbm.simplify(brt.m)
lapply(brt.m.simplify$pred.list,function(x) vars[x])

brt.m.predictors <- vars[brt.m.simplify$pred.list[[7]]]
brt.m.formula    <- log(m)~family+log(k)+log(amax)

brt.m <- Brter(brt.m.formula,exp,tree.complexity = 10,learning.rate = 0.002,bag.fraction = 0.5,ntrees = 3050)
#brt.m$fit(fb)

brt.m.cv <- brt.m$cross(fb,folds=100)

dfr <- cbind(par='m',brt.m.cv$folds[,c('hat','obs')])

ggplot(dfr) + 
  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  theme_bw(base_size=20) + 
  labs(x='Observed',y='Predicted',title=as.character(brt.m.formula)[3])

saver(brt.m,brt.m.formula,brt.m.predictors,brt.m.cv,name='brt_m')

# let's try using the pars argument
#
#brt.test <- Brter(formula,exp,ntrees=5000,learning.rate=0.001)
#
#brt.test$fit(fb)
#pars <- brt.test$brt$gbm.call$gbm.x
#
#res <- list()
#pars.update <- pars
#
#for(i in 1:(length(pars)-1)) {
#  
#  brt.test$fit(fb,pars.update)
#  cat('predictors:',brt.test$brt$gbm.call$predictor.names,'\n')
#  
#  rinfl <- relative.influence(brt.test$brt,n.trees = brt.test$brt$gbm.call$best.trees)
#  
#  res[[i]] <- list()
#  res[[i]][['summary']]   <- brt.test$cross(fb,pars=pars.update)$summary 
#  res[[i]][['influence']] <- data.frame(predictor=brt.test$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
#  res[[i]][['drop']]      <- brt.test$brt$gbm.call$predictor.names[which.min(rinfl)]
#  
#  pars.update <- pars.update[-which.min(rinfl)]
#  
#}
#
# results are consistent with gbm.simplify()
#unlist(lapply(res,function(x) x$drop))
#

par.names <- c("class","order","family","sex","trophic","habit","log(temp)","log(lmax)","log(k)","log(amax)")
pars.update <- 1:length(par.names)
  
brt.m.res <- list()

for(i in 1:(length(par.names)-1)) {
  
 formula.update <- as.formula(paste('log(m)~',paste(par.names[pars.update],collapse='+'),sep=''))
 
 brt.update <- Brter(formula.update,exp,ntrees=5000,learning.rate=0.001)
 
 brt.update$fit(fb)
 
 rinfl <- relative.influence(brt.update$brt,n.trees = brt.update$brt$gbm.call$best.trees)
 cross <- brt.update$cross(fb,folds=100) 
 
 brt.m.res[[i]] <- list()
 brt.m.res[[1]][['formula']]    <- formula.update
 brt.m.res[[i]][['summary']]    <- cross$summary 
 brt.m.res[[i]][['folds']]      <- cross$folds
 brt.m.res[[i]][['influence']]  <- data.frame(predictor=brt.update$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
 brt.m.res[[i]][['drop']]       <- brt.update$brt$gbm.call$predictor.names[which.min(rinfl)]
 
 pars.update <- pars.update[-which.min(rinfl)]
 
}

saver(brt.m.res,name='brt_m_res')

# FIGURES
par.drop <- unlist(lapply(brt.m.res,function(x) x$drop))

dfr <- data.frame()
dfr <- rbind(dfr,data.frame(id=1,order=1:9,predictor.removed=c('none',par.drop[-length(par.drop)]),value=unlist(lapply(brt.m.res,function(x) {y<-x$summary['mpe','mean'];y})),label='Mean prediction error'))
dfr <- rbind(dfr,data.frame(id=1,order=1:9,predictor.removed=c('none',par.drop[-length(par.drop)]),value=unlist(lapply(brt.m.res,function(x) {y<-x$summary['dev','mean'];y})),label='Mean prediction deviance'))

ggplot(dfr) + 
  geom_line(aes(x=order,y=value,group=id),size=2) + 
  scale_x_discrete(labels=dfr$predictor.removed) + 
  facet_wrap(~label,ncol=1,scale='free_y') +
  labs(x='Predictor removed',y='') +
  theme_bw(base_size=20) 

   