
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


############################
# fishbase                 #
############################

#formula <- log(m)~family+sex+trophic+habit+log(temp)+log(lmat)+log(linf)+log(k)+log(amax)+log(amat)
#vars    <- all.vars(formula)
#frame <- model.frame(formula,fb,na.action=na.pass)
#names(frame) <- c("log.m","family","sex","trophic","habit","log.temp","log.lmat","log.linf","log.k","log.amax","log.amat")
#brt.m <- gbm.step(data = frame[!is.na(frame$log.m),],
#                  gbm.y = 1,
#                  gbm.x = 2:ncol(frame), 
#                  family = "gaussian",
#                  tree.complexity = 10,
#                  learning.rate = 0.001,
#                  bag.fraction = 0.5,
#                  max.trees = 10000)
#
#summary(brt.m)
#
#brt.m.simplify <- gbm.simplify(brt.m)
#lapply(brt.m.simplify$pred.list,function(x) vars[x])
#
#brt.m.predictors <- vars[brt.m.simplify$pred.list[[7]]]
#brt.m.formula    <- log(m)~family+log(k)+log(amax)
#
#brt.m <- Brter(brt.m.formula,exp,tree.complexity = 10,learning.rate = 0.002,bag.fraction = 0.5,ntrees = 3050)
#brt.m$fit(fb)
#
#brt.m.cv <- brt.m$cross(fb,folds=100)
#
#dfr <- cbind(par='m',brt.m.cv$folds[,c('hat','obs')])
#
#ggplot(dfr) + 
#  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
#  geom_abline(a=0,b=1) + 
#  theme_bw(base_size=20) + 
#  labs(x='Observed',y='Predicted',title=as.character(brt.m.formula)[3])
#
#saver(brt.m,brt.m.formula,brt.m.predictors,brt.m.cv,name='brt_m')

# preliminary screening removed class, order and amat
# inclusion of lmat appeared to inflate the minimum deviance and prediction error
par.names <- c("family","trophic","habit","sex","log(temp)","log(linf)","log(k)","log(amax)")
pars.update <- 1:length(par.names)
  
res <- list()

for(i in 1:(length(par.names)-1)) {
  
 formula.update <- as.formula(paste('log(m)~',paste(par.names[pars.update],collapse='+'),sep=''))
 
 brt.update <- Brter(formula.update,exp,ntrees=5000,learning.rate=0.001)
 
 brt.update$fit(fb)
 
 rinfl <- relative.influence(brt.update$brt,n.trees = brt.update$brt$gbm.call$best.trees)
 cross <- brt.update$cross(fb,folds=10) 
 
 res[[i]] <- list()
 res[[i]][['formula']]    <- formula.update
 res[[i]][['summary']]    <- cross$summary 
 res[[i]][['folds']]      <- cross$folds
 res[[i]][['influence']]  <- data.frame(predictor=brt.update$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
 res[[i]][['drop']]       <- brt.update$brt$gbm.call$predictor.names[which.min(rinfl)]
 
 pars.update <- pars.update[-which.min(rinfl)]
 
}

# FIGURES
par.drop <- unlist(lapply(res,function(x) x$drop))
par.remove <- c('none',par.drop[-length(par.drop)])

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
 dfr1 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=res[[i]]$folds[,'mpe'],label='Mean prediction error')
 dfr2 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=res[[i]]$folds[,'dev'],label='Mean prediction deviance')
 dfr  <- rbind(dfr,dfr1,dfr2)
}

fig <- ggplot(dfr,aes(x=order,y=value,group=id)) +
  stat_summary(fun.ymin=function(x) mean(x) - sd(x)/sqrt(10),fun.ymax=function(x) mean(x) + sd(x)/sqrt(10),geom='ribbon',alpha=0.2) +
  stat_summary(fun.y=function(x) mean(x),geom='line',lwd=1.5) +
  scale_x_discrete(labels=par.remove) + 
  facet_wrap(~label,ncol=1,scale='free_y') +
  labs(x='Predictor removed',y='') +
  theme_bw(base_size=20)

pdfr(fig,width=12,name='fb_v2/m_brt_res_fig1')

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
  dfr  <- rbind(dfr,data.frame(id=1,order=i,fold=1:10,predictor.removed=paste('drop ',i-1,': ',par.remove[i],sep=''),obs=res[[i]]$folds[,'obs'],hat=res[[i]]$folds[,'hat']))
}


fig <- ggplot(dfr,aes(x=obs,y=hat)) +
  geom_point(size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~predictor.removed) + 
  theme_bw(base_size=20) +
  labs(x='Observed value',y='Prediction')

pdfr(fig,name='fb_v2/m_brt_res_fig2')


# SAVE RESULTS
brt.initial <- Brter(log(m)~family+trophic+habit+sex+log(temp)+log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)
brt.final   <- Brter(log(m)~family+log(temp)+log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)

brt.initial$fit(fb)
brt.final$fit(fb)

#summary(brt.final$brt)

saver(res,brt.initial,brt.final,name='fb_v2/m_brt_res')


############################
# gislason                 #
############################

# preliminary screening removed class, order and amat
# inclusion of lmat appeared to inflate the minimum deviance and prediction error
par.names <- c("family","trophic","habit","sex","log(temp)","log(linf)","log(k)","log(amax)")
pars.update <- 1:length(par.names)

res <- list()

for(i in 1:(length(par.names)-1)) {
  
  formula.update <- as.formula(paste('log(m)~',paste(par.names[pars.update],collapse='+'),sep=''))
  
  brt.update <- Brter(formula.update,exp,ntrees=5000,learning.rate=0.001)
  
  brt.update$fit(gs)
  
  rinfl <- relative.influence(brt.update$brt,n.trees = brt.update$brt$gbm.call$best.trees)
  cross <- brt.update$cross(fb,folds=10) 
  
  res[[i]] <- list()
  res[[i]][['formula']]    <- formula.update
  res[[i]][['summary']]    <- cross$summary 
  res[[i]][['folds']]      <- cross$folds
  res[[i]][['influence']]  <- data.frame(predictor=brt.update$brt$gbm.call$predictor.names,influence=as.numeric(rinfl/sum(rinfl)))
  res[[i]][['drop']]       <- brt.update$brt$gbm.call$predictor.names[which.min(rinfl)]
  
  pars.update <- pars.update[-which.min(rinfl)]
  
}

# FIGURES
par.drop <- unlist(lapply(res,function(x) x$drop))
par.remove <- c('none',par.drop[-length(par.drop)])

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
  dfr1 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=res[[i]]$folds[,'mpe'],label='Mean prediction error')
  dfr2 <- data.frame(id=1,order=i,fold=1:10,predictor.removed=par.remove[i],value=res[[i]]$folds[,'dev'],label='Mean prediction deviance')
  dfr  <- rbind(dfr,dfr1,dfr2)
}

fig <- ggplot(dfr,aes(x=order,y=value,group=id)) +
  stat_summary(fun.ymin=function(x) mean(x) - sd(x)/sqrt(10),fun.ymax=function(x) mean(x) + sd(x)/sqrt(10),geom='ribbon',alpha=0.2) +
  stat_summary(fun.y=function(x) mean(x),geom='line',lwd=1.5) +
  scale_x_discrete(labels=par.remove) + 
  facet_wrap(~label,ncol=1,scale='free_y') +
  labs(x='Predictor removed',y='') +
  theme_bw(base_size=20)

pdfr(fig,width=12,name='gs_v2/m_brt_res_fig1')

dfr <- data.frame()
for(i in 1:(length(par.names)-1)) {
  
  dfr  <- rbind(dfr,data.frame(id=1,order=i,fold=1:10,predictor.removed=paste('drop ',i-1,': ',par.remove[i],sep=''),obs=res[[i]]$folds[,'obs'],hat=res[[i]]$folds[,'hat']))
}


fig <- ggplot(dfr,aes(x=obs,y=hat)) +
  geom_point(size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~predictor.removed) + 
  theme_bw(base_size=20) +
  labs(x='Observed value',y='Prediction')

pdfr(fig,name='gs_v2/m_brt_res_fig2')


# FINAL MODEL
brt.initial <- Brter(log(m)~family+trophic+habit+sex+log(temp)+log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)
brt.final   <- Brter(log(m)~family+log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)

brt.initial$fit(gs)
brt.final$fit(gs)

#summary(brt.final$brt)

saver(res,brt.initial,brt.final,name='gs_v2/m_brt_res')

#################
# SUMMARY TABLE #
#################

rm(list=ls())
source('collate.R')
source('load_data.R')
source('utils.R')

dfr <- data.frame()

loader('fb_v2/m_brt_res')
brt.cv <- brt.initial$cross(fb)
dfr <- rbind(dfr,
             data.frame(source='(7a) initial',db='fb',n=brt.initial$n(fb),mpe=round(brt.cv$summary['mpe',1],2),dev=round(brt.cv$summary['dev',1],2))
)
brt.cv <- brt.final$cross(fb)
dfr <- rbind(dfr,
             data.frame(source='(7a) final',db='fb',n=brt.final$n(fb),mpe=round(brt.cv$summary['mpe',1],2),dev=round(brt.cv$summary['dev',1],2))
)
windows(width=4)
summary(brt.final$brt)
savePlot(file='C:/PROJECTS/FISHNETS/res/fb_v2/m_brt_res_fig3.pdf',type='pdf')
dev.off()

loader('gs_v2/m_brt_res')
brt.cv <- brt.initial$cross(gs)
dfr <- rbind(dfr,
             data.frame(source='(7b) initial',db='gs',n=brt.initial$n(gs),mpe=round(brt.cv$summary['mpe',1],2),dev=round(brt.cv$summary['dev',1],2))
)
brt.cv <- brt.final$cross(gs)
dfr <- rbind(dfr,
             data.frame(source='(7b) final',db='gs',n=brt.final$n(gs),mpe=round(brt.cv$summary['mpe',1],2),dev=round(brt.cv$summary['dev',1],2))
)
windows(width=4)
summary(brt.final$brt)
savePlot(file='C:/PROJECTS/FISHNETS/res/gs_v2/m_brt_res_fig3.pdf',type='pdf')
dev.off()


write.csv(dfr,file='C:/PROJECTS/FISHNETS/res/cvbrt_v2.csv')


