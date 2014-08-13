
############################
# boosted regression trees #
############################

library(ggplot2)

rm(list=ls())
source('collate.R')
source('load_data.R')
source('utils.R')

# according to life history theory, primary correlates
# with m are: k and amax; via the BH invariants
# m/k and m*amax. Therefore need to optimise estimation
# of these too.

##############
# groom data #
##############

fb[which(fb$temp<=0),'temp'] <- NA
fb[which(fb$m>2),'m']        <- NA

gs[which(gs$m>2),'m']        <- NA

# calculate amat in fishbase using VB growth equation
obj <- function(alpha) lmat - linf * (1 - exp(-k * (alpha - t0)))

for(i in 1:length(fb$amat)) {
  
  if(is.na(fb$amat[i])) {
    
    lmat <- fb$lmat[i]
    linf <- fb$linf[i]
    k    <- fb$k[i]
    t0   <- fb$t0[i]
    
    if(any(is.na(c(lmat,linf,k,t0)))) next
    if(lmat>linf) next
    if(lmat<(linf * (1 - exp(-k * (0 - t0))))) next
    
    fb$amat[i] <- uniroot(obj,interval=c(0,100))$root
    
  }
}

# estimate amat in gislasson database using VB growth
# equation and assuming t0 = 0 (Gislason 2010)
gs$amat <- NA

for(i in 1:length(gs$amat)) {
  
  if(is.na(gs$amat[i])) {
    
    lmat <- gs$lmat[i]
    linf <- gs$linf[i]
    k    <- gs$k[i]
    t0   <- 0
    
    if(any(is.na(c(lmat,linf,k,t0)))) next
    if(lmat>linf) next
    if(lmat<(linf * (1 - exp(-k * (0 - t0))))) next
    
    gs$amat[i] <- uniroot(obj,interval=c(0,100))$root
    
  }
}


############################
# fishbase                 #
############################

formula <- log(m)~family+sex+trophic+habit+log(temp)+log(lmat)+log(linf)+log(k)+log(amax)+log(amat)
vars    <- all.vars(formula)
frame <- model.frame(formula,fb,na.action=na.pass)
names(frame) <- c("log.m","family","sex","trophic","habit","log.temp","log.lmat","log.linf","log.k","log.amax","log.amat")
brt.m <- gbm.step(data = frame[!is.na(frame$log.m),],
                  gbm.y = 1,
                  gbm.x = 2:ncol(frame), 
                  family = "gaussian",
                  tree.complexity = 10,
                  learning.rate = 0.001,
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
#brt.test <- Brter(log(m)~family+trophic+habit+sex+log(temp)+log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)
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


# preliminary screening removed class, order and amat
# inclusion of lmat appeared to inflate the minimum deviance and prediction error
par.names <- c("family","trophic","habit","sex","log(temp)","log(linf)","log(k)","log(amax)")
pars.update <- 1:length(par.names)
  
brt.m.res <- list()

for(i in 1:(length(par.names)-1)) {
  
 formula.update <- as.formula(paste('log(m)~',paste(par.names[pars.update],collapse='+'),sep=''))
 
 brt.update <- Brter(formula.update,exp,ntrees=5000,learning.rate=0.001)
 
 brt.update$fit(fb)
 
 rinfl <- relative.influence(brt.update$brt,n.trees = brt.update$brt$gbm.call$best.trees)
 cross <- brt.update$cross(fb,folds=10) 
 
 brt.m.res[[i]] <- list()
 brt.m.res[[i]][['formula']]    <- formula.update
 brt.m.res[[i]][['summary']]    <- cross$summary 
 brt.m.res[[i]][['folds']]      <- cross$folds
 brt.m.res[[i]][['influence']]  <- data.frame(predictor=brt.update$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
 brt.m.res[[i]][['drop']]       <- brt.update$brt$gbm.call$predictor.names[which.min(rinfl)]
 
 pars.update <- pars.update[-which.min(rinfl)]
 
}

saver(brt.m.res,name='fb/m_brt_res')

# FIGURES
par.drop <- unlist(lapply(brt.m.res,function(x) x$drop))
par.remove <- c('none',par.drop[-length(par.drop)])

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
 dfr1 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=brt.m.res[[i]]$folds[,'mpe'],label='Mean prediction error')
 dfr2 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=brt.m.res[[i]]$folds[,'dev'],label='Mean prediction deviance')
 dfr  <- rbind(dfr,dfr1,dfr2)
}

fig <- ggplot(dfr,aes(x=order,y=value,group=id)) +
  stat_summary(fun.ymin=function(x) mean(x) - sd(x)/sqrt(10),fun.ymax=function(x) mean(x) + sd(x)/sqrt(10),geom='ribbon',alpha=0.2) +
  stat_summary(fun.y=function(x) mean(x),geom='line',lwd=1.5) +
  scale_x_discrete(labels=par.remove) + 
  facet_wrap(~label,ncol=1,scale='free_y') +
  labs(x='Predictor removed',y='') +
  theme_bw(base_size=20)

pdfr(fig,width=12,name='fb/m_brt_res_fig1')

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
  dfr  <- rbind(dfr,data.frame(id=1,order=i,fold=1:10,predictor.removed=paste('drop ',i-1,': ',par.remove[i],sep=''),obs=brt.m.res[[i]]$folds[,'obs'],hat=brt.m.res[[i]]$folds[,'hat']))
}


fig <- ggplot(dfr,aes(x=obs,y=hat)) +
  geom_point(size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~predictor.removed) + 
  theme_bw(base_size=20) +
  labs(x='Observed value',y='Prediction')

pdfr(fig,name='fb/m_brt_res_fig2')


# FINAL MODEL
formula.final <-log(m)~family+log(temp)+log(linf)+log(k)+log(amax)

brt.final <- Brter(formula.final,exp,ntrees=5000,learning.rate=0.001)

brt.final$fit(fb)

summary(brt.final$brt)

# performance statistics
brt.m.res[[4]][['formula']]
brt.m.res[[i]][['summary']]


############################
# gislason                 #
############################

# preliminary screening removed class, order and amat
# inclusion of lmat appeared to inflate the minimum deviance and prediction error
par.names <- c("family","trophic","habit","sex","log(temp)","log(linf)","log(k)","log(amax)")
pars.update <- 1:length(par.names)

brt.m.res <- list()

for(i in 1:(length(par.names)-1)) {
  
  formula.update <- as.formula(paste('log(m)~',paste(par.names[pars.update],collapse='+'),sep=''))
  
  brt.update <- Brter(formula.update,exp,ntrees=5000,learning.rate=0.001)
  
  brt.update$fit(gs)
  
  rinfl <- relative.influence(brt.update$brt,n.trees = brt.update$brt$gbm.call$best.trees)
  cross <- brt.update$cross(fb,folds=10) 
  
  brt.m.res[[i]] <- list()
  brt.m.res[[i]][['formula']]    <- formula.update
  brt.m.res[[i]][['summary']]    <- cross$summary 
  brt.m.res[[i]][['folds']]      <- cross$folds
  brt.m.res[[i]][['influence']]  <- data.frame(predictor=brt.update$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
  brt.m.res[[i]][['drop']]       <- brt.update$brt$gbm.call$predictor.names[which.min(rinfl)]
  
  pars.update <- pars.update[-which.min(rinfl)]
  
}

saver(brt.m.res,name='gs/m_brt_res')

# FIGURES
par.drop <- unlist(lapply(brt.m.res,function(x) x$drop))
par.remove <- c('none',par.drop[-length(par.drop)])

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
  dfr1 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=brt.m.res[[i]]$folds[,'mpe'],label='Mean prediction error')
  dfr2 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=brt.m.res[[i]]$folds[,'dev'],label='Mean prediction deviance')
  dfr  <- rbind(dfr,dfr1,dfr2)
}

fig <- ggplot(dfr,aes(x=order,y=value,group=id)) +
  stat_summary(fun.ymin=function(x) mean(x) - sd(x)/sqrt(10),fun.ymax=function(x) mean(x) + sd(x)/sqrt(10),geom='ribbon',alpha=0.2) +
  stat_summary(fun.y=function(x) mean(x),geom='line',lwd=1.5) +
  scale_x_discrete(labels=par.remove) + 
  facet_wrap(~label,ncol=1,scale='free_y') +
  labs(x='Predictor removed',y='') +
  theme_bw(base_size=20)

pdfr(fig,width=12,name='gs/m_brt_res_fig1')

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
  dfr  <- rbind(dfr,data.frame(id=1,order=i,fold=1:10,predictor.removed=paste('drop ',i-1,': ',par.remove[i],sep=''),obs=brt.m.res[[i]]$folds[,'obs'],hat=brt.m.res[[i]]$folds[,'hat']))
}


fig <- ggplot(dfr,aes(x=obs,y=hat)) +
  geom_point(size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~predictor.removed) + 
  theme_bw(base_size=20) +
  labs(x='Observed value',y='Prediction')

pdfr(fig,name='gs/m_brt_res_fig2')


# FINAL MODEL
formula.final <-log(m)~family+log(linf)+log(k)+log(amax)

brt.final <- Brter(formula.final,exp,ntrees=5000,learning.rate=0.001)

brt.final$fit(gs)

summary(brt.final$brt)

   