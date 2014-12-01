
########################
# single-node fishnets #
########################

library(ggplot2)

rm(list=ls())
source('collate.R')
source('load_data.R')
source('utils.R')
source('cvfishnet.R')



#########################
# Node 1 - Charnov 1990 #
#########################

mfishnet <- Fishnet(

  m         = MCharnov1990Fitted()
  
)

# cross validate fishnet

# fishbase
data.fb <- model.frame('m~species+genus+family+order+class+amat',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet1_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet1_prediction')

# gislason
data.gs <- model.frame('m~species+genus+family+order+class+amat',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet1_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet1_prediction')


#########################
# Node 2 - Charnov 1993 #
#########################

mfishnet <- Fishnet(
    
  m         = MCharnov1993Fitted()
  
)

# cross validate fishnet

# fishbase
data.fb <- model.frame('m~species+genus+family+order+class+k',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet2_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet2_prediction')

# gislason
data.gs <- model.frame('m~species+genus+family+order+class+k',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet2_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet2_prediction')

#########################
# Node 3 - Hoenig 1983 #
#########################

mfishnet <- Fishnet(
    
  m         = MHoenig1983Fitted()
  
)

# cross validate fishnet

# fishbase
data.fb <- model.frame('m~species+genus+family+order+class+amax',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet3_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet3_prediction')

# gislason
data.gs <- model.frame('m~species+genus+family+order+class+amax',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet3_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet3_prediction')

#########################
# Node 4 - Quinn 1999   #
#########################

mfishnet <- Fishnet(
    
  m         = MQuinn1999Fitted()
  
)

# cross validate fishnet

# fishbase
data.fb <- model.frame('m~species+genus+family+order+class+amax',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet4_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet4_prediction')

# gislason
data.gs <- model.frame('m~species+genus+family+order+class+amax',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet4_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet4_prediction')


#########################
# Node 5 - Charnov 2013 #
#########################

mfishnet <- Fishnet(
    
  m         = MCharnov2013Fitted()
  
)

# cross validate fishnet

# fishbase
data.fb <- model.frame('m~species+genus+family+order+class+lmat+linf+k',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet5_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet5_prediction')

# gislason
data.gs <- model.frame('m~species+genus+family+order+class+lmat+linf+k',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet5_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet5_prediction')


#########################
# Node 6 - Pauly 1980   #
#########################

mfishnet <- Fishnet(
    
  m         = MPauly1980Fitted()
  
)

# cross validate fishnet

# fishbase
data.fb <- model.frame('m~species+genus+family+order+class+temp+linf+k',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet6_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet6_prediction')

# gislason
data.gs <- model.frame('m~species+genus+family+order+class+temp+linf+k',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet6_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet6_prediction')


#########################
# Node 7a - BRT fb      #
#########################

mfishnet <- Fishnet(

  m         = Brter(log(m)~log(temp)+log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)
  
)

# cross validate fishnet

data.fb <- model.frame('m~species+genus+family+order+class+temp+linf+k+amax',fb)
data.fb <- data.fb[!is.na(data.fb$m),]

res <- cvfishnet(mfishnet,data.fb,byspecies=T)
saver(mfishnet,res,name='fb_v2/mfishnet7a_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n')+
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='fb_v2/mfishnet7a_prediction')


#########################
# Node 7b - BRT gs      #
#########################

mfishnet <- Fishnet(
  
  m         = Brter(log(m)~log(linf)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)
  
)

# cross validate fishnet

data.gs <- model.frame('m~species+genus+family+order+class+linf+k+amax',gs)
data.gs <- data.gs[!is.na(data.gs$m),]

res <- cvfishnet(mfishnet,data.gs,byspecies=T)
saver(mfishnet,res,name='gs_v2/mfishnet7b_res')

dfr <- data.frame(res$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

fig <- ggplot(subset(dfr,label=='Predicted Values'),aes(x=obs,y=hat)) + 
    geom_point(size=4,alpha=0.3) + 
    geom_abline(aes(intercept=0,slope=1),colour="#990000", linetype="dashed",size=1) +
    theme_bw(base_size=25) +
    labs(x='Observed',y='Predicted\n') +
    scale_x_continuous(limits=c(0, 2)) +
    scale_y_continuous(limits=c(0, 2))

pdfr(fig,width=10,name='gs_v2/mfishnet7b_prediction')


#################
# SUMMARY TABLE #
#################

rm(list=ls())
source('collate.R')
source('load_data.R')
source('utils.R')

dfr <- data.frame()

loader('fb_v2/mfishnet1_res')
dfr <- rbind(dfr,
  data.frame(source='(1)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)
loader('gs_v2/mfishnet1_res')
dfr <- rbind(dfr,
  data.frame(source='(1)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('fb_v2/mfishnet2_res')
dfr <- rbind(dfr,
             data.frame(source='(2)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)
loader('gs_v2/mfishnet2_res')
dfr <- rbind(dfr,
             data.frame(source='(2)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('fb_v2/mfishnet3_res')
dfr <- rbind(dfr,
             data.frame(source='(3)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)
loader('gs_v2/mfishnet3_res')
dfr <- rbind(dfr,
             data.frame(source='(3)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('fb_v2/mfishnet4_res')
dfr <- rbind(dfr,
             data.frame(source='(4)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)
loader('gs_v2/mfishnet4_res')
dfr <- rbind(dfr,
             data.frame(source='(4)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('fb_v2/mfishnet5_res')
dfr <- rbind(dfr,
             data.frame(source='(5)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)
loader('gs_v2/mfishnet5_res')
dfr <- rbind(dfr,
             data.frame(source='(5)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('fb_v2/mfishnet6_res')
dfr <- rbind(dfr,
             data.frame(source='(6)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)
loader('gs_v2/mfishnet6_res')
dfr <- rbind(dfr,
             data.frame(source='(6)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('fb_v2/mfishnet7a_res')
dfr <- rbind(dfr,
             data.frame(source='(7a)',db='fb',n=mfishnet$nodes[['m']]$n(fb),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

loader('gs_v2/mfishnet7b_res')
dfr <- rbind(dfr,
             data.frame(source='(7b)',db='gs',n=mfishnet$nodes[['m']]$n(gs),mpe=round(res$summary['mpe',1],2),dev=round(res$summary['dev',1],2))
)

write.csv(dfr,file='C:/PROJECTS/FISHNETS/res/cvfishnet_v2.csv')


