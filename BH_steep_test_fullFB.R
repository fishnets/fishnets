# steepness tests;

# Source in the package ----
source('collate.R')

# Load the Fishbase data
fb <- FishbaseWeb$read('data/fishbase-web')
# Limit to the 7940 with both k and linf
fb <- subset(fb,!is.na(k) & !is.na(linf))
# An an id column for indexing later
fb$id <- 1:nrow(fb)
# Add a dummy row for helping with predictor nodes
# that need to have at least two predictors
fb$dummy <- 1.0


steep <- Steepness$create('./data/steepness')
steep_merged <- Steepness$merge(steep,fb)

# groom
# fecundity of 1 is absurd unless boolean
steep_merged$fecundity[steep_merged$fecundity==1] <- NA
steep_merged$recsigma <- NA
# build a net for steepnes. Use Bayesian nodes in an attempt to not overfit (i.e., to egt better predictive power)

logit <- function(x) log(x/(1-x))
logit_inv <- function(xt) 1/(1 + exp(-xt))
BH_tr <- function(h) h/0.8-0.25
BH_tr_inv <- function(ht) (ht+0.25)*0.8
logit_BH <- function(h) logit(BH_tr(h))
logit_BH_inv <- function(ht) BH_tr_inv(logit_inv(ht))

BH_net_full <- Fishnet(
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit'),
  depthmax  = TaxonomicImputer('depthmax',c(log,exp),5),
  trophic   = TaxonomicImputer('trophic',c(log,exp),3),
  lmax      = TaxonomicImputer('lmax',c(log,exp),5),
  amax      = TaxonomicImputer('amax',c(log,exp),5),
  
  linf      = Bayser(log(linf) ~ f(family,model="iid")+f(class,model="iid")+log(lmax),exp),
  fecundity = Bayser(log(fecundity) ~ f(family,model="iid")+f(class,model="iid") + log(linf) + log(depthmax),exp),
  k         = Bayser(log(k) ~ f(family,model="iid") + log(linf) + f(habit,model="iid") + log(depthmax),exp),
  m         = Bayser(log(m) ~ f(family,model="iid")+f(class,model="iid")+log(k)+log(linf)+f(habit,model="iid")+log(depthmax)+trophic,exp),
  lmat      = Bayser(log(lmat) ~ f(family,model="iid")+log(k)+log(linf)+f(habit,model="iid")+log(depthmax),exp),
  recsigma  = RecsigmaThorsonEtAl2014(),
  mean_BH_z  = Brter(logit_BH(mean_BH_z) ~  habit + log(linf) + log(k) + log(m)+ log(fecundity) +recsigma + trophic+log(depthmax),transform = logit_BH_inv,ntrees =0,bag.fraction=0.9)
  
)

# fit the BH_et to the summarised fishbase data
BH_net_full$fit(steep_merged,impute = T)

# function to make testset for cross validation
make_testset <- function(net,data,name){
  testset <- data.frame(net$data[,-which(colnames(net$data) == name)], name = data[name])
  testset
}

testset <- make_testset(BH_net_full,steep_merged,'mean_BH_z')

# CV meta function
crossval <- function(net,data,name){
  net$nodes[[name]]$cross(make_testset(net,data,name))
}

# do CV on nodes other than steepness
crossval(BH_net_full,steep_red,'linf')
crossval(BH_net_full,steep_red,'lmat')


# # data with only imputed steepness values, no originals
# trainset <- BH_net_full$data
# trainset$mean_BH_z[!is.na(steep_red$mean_BH_z)] <- NA

BH_net_full$nodes$mean_BH_z  = Brter(logit_BH(mean_BH_z) ~ log(linf) + log(k) + log(m)+ log(fecundity) +recsigma + trophic+log(depthmax),transform = logit_BH_inv,ntrees =3000,bag.fraction=0.8)


BH_net_full$nodes$mean_BH_z$fit(testset)

plot(BH_net_full$nodes$mean_BH_z$brt)
summary(BH_net_full$nodes$mean_BH_z$brt)

Predicted <- BH_net_full$nodes$mean_BH_z$predict(testset)[!is.na(steep_merged$mean_BH_z)]
Observed <- steep_merged$mean_BH_z[!is.na(steep_merged$mean_BH_z)]

self_pred <- lm(logit_BH(Observed)~logit_BH(Predicted))
summary(self_pred)

plot(logit_BH(Predicted),logit_BH(Observed),pch=16)
abline(self_pred$coeff[1],self_pred$coeff[2],col=2,lwd=2)
abline(a=0,b=1,lwd=2)

# it's cheating a little
# jacknifing - could be a node feature in fishnets alongside CV

jacknife_cv <- function(data,net,node){
  testnet <- net
  data = data[!is.na(data[[node]]),]
  
  pred <- vector(,nrow(data))
  for (i in 1:nrow(data)){
    cat('CV for observation ',i,'\n')
    train <- data[-i,]
    test <- data[i,]
    test[[node]] <- NA
    testnet$nodes[[node]]$fit(train)
    pred[i] <- testnet$nodes[[node]]$predict(test)
    
  }
  data.frame(Predicted = pred,Observed = data[[node]])
}

steep_cv <- jacknife_cv(testset,BH_net_full,'mean_BH_z')

lm_pred_steep <- lm(Observed~Predicted,data=steep_cv)
summary(lm_pred_steep)

plot(steep_cv,pch=16)
abline(lm_pred_steep$coeff[1],lm_pred_steep$coeff[2],col=2,lwd=2)
abline(0,1,lwd=2)

BH_net_full$nodes$mean_BH_z  = Bayser(logit_BH(mean_BH_z) ~ f(species,model='iid') +log(linf) + log(k) + log(m)+ log(fecundity) +log(recsigma)+log(m)*log(recsigma) + log(trophic) + log(depthmax),logit_BH_inv)

BH_net_full$nodes$mean_BH_z$fit(testset)

steep_cv_bayes <- jacknife_cv(testset,BH_net_full,'mean_BH_z')

lm_pred_steep <- lm(Observed~Predicted,data=steep_cv_bayes)
summary(lm_pred_steep)

plot(steep_cv_bayes,pch=16)
abline(lm_pred_steep$coeff[1],lm_pred_steep$coeff[2],col=2,lwd=2)
abline(0,1,lwd=2)


####### Bluenose --------

bwa <- BH_net_full$sample(list(
  species = 'Hyperoglyphe antarctica',
  # Maximum length, temperature and 
  # maximum depth from Fishbase
  lmax = 140,
  temp = 11,
  depthmax = 1500,
  # Female growth and max age from 
  # Horn et al 2010
  linf = 92.5,
  k = 0.071,
  amax = 71  
),samples = 1000)

ggplot(bwa) + 
  geom_bar(aes(x=mean_BH_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0.2,1)) + 
  labs(x='Steepness (z)',y='Density')

ggplot(bwa) + 
  geom_point(aes(x=m,y=mean_BH_z),alpha=0.4) + 
  scale_x_log10(breaks=seq(0.1,1.1,0.2)) + 
  scale_y_log10(breaks=seq(0.1,1.1,0.2)) + 
  labs(x='Natural mortality',y='Steepness')

# how much information is gained from life-history

bwa.org <- BH_net_full$sample(list(
  species = 'Hyperoglyphe antarctica'),samples = 10000)

ggplot(bwa.org) + 
  geom_bar(aes(x=mean_BH_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0.2,1)) + 
  labs(x='Steepness (z)',y='Density')


# Comparison of estiamted vs 'data' -----

BH_net_full_test <- BH_net_full

#' Plot density histograms
plot_samples <- function(samples,inp_data,species_,pars=c('linf','k','m','mean_BH_z')){
  datas = subset(inp_data,species==species_)
  
  melted <- melt(samples[,pars])
  data_melted <- melt(datas[,pars])
  ggplot(melted,aes(x=value)) +
    geom_bar(data=data_melted,aes(y = 1)) +
    geom_density(fill=hsv(0,0.7,0.7),alpha=0.5) +
    facet_wrap(~variable,scales='free') + 
    labs(x='',y='Density') + 
    theme(strip.text.x=element_text(size=10))
}


# fit test net
BH_net_full_test$fit(subset(steep_red,species!='Gadus morhua'),impute = T)

# predictions

preds.nlh <- BH_net_full_test$sample(list(
  species = 'Gadus morhua'
),samples=1000)

plot_samples(preds.nlh,steep_merged,'Gadus morhua')

preds.slh <- BH_net_full_test$sample(list(
  species = 'Gadus morhua',
  swimmode = 'subcarangiform',
  habit = 'benthopelagic',
  depthmax = 600,
  temp = 5,
  lmax = 132
),samples = 1000)

plot_samples(preds.slh,steep_merged,
  'Gadus morhua'
)

preds.lh <- BH_net_full$sample(dists(
  species =  Fixed('Gadus morhua'),
  swimmode = Fixed('subcarangiform'),
  habit = Fixed('benthopelagic'),
  depthmax = Fixed(600),
  temp = Trapezoid(2,3,6,7),
  lmax = Fixed(132),
  linf = Normal(110,20),
  k = Triangle(0.07,0.13,0.35),
  amax=Fixed(20)
),1000)


plot_samples(preds.lh,steep_merged,'Gadus morhua')

# without imputation

BH_net_full_noImputation <- BH_net_full_test
BH_net_full_noImputation$fit(subset(steep_merged,species!='Gadus morhua'), impute = F)

preds_noImputation <- BH_net_full_noImputation$sample(dists(
  species =  Fixed('Gadus morhua'),
  swimmode = Fixed('subcarangiform'),
  habit = Fixed('benthopelagic'),
  depthmax = Fixed(600),
  temp = Trapezoid(2,3,6,7),
  lmax = Fixed(132),
  linf = Normal(110,20),
  k = Triangle(0.07,0.13,0.35),
  amax=Fixed(20)
),1000)

plot_samples(preds_noImputation,'Gadus morhua')

###### Skipjack ------

BH_net_full_test$fit(subset(steep_merged,species!='Katsuwonus pelamis'))

plot_samples(
  
  BH_net_full$sample(list(
    species = 'Katsuwonus pelamis',
    family = 'Scombridae'
  )),
  
  'Katsuwonus pelamis'
)

plot_samples(
  
  BH_net_full$sample(dists(
    species = Fixed('Katsuwonus pelamis'),
    family = Fixed('Scombridae'),
    depthmax = Fixed(260),
    lmax = Fixed(90.5)
  )),
  
  'Katsuwonus pelamis'
)


plot_samples(
  
  BH_net_full$sample(dists(
    species = Fixed('Katsuwonus pelamis'),
    family = Fixed('Scombridae'),
    depthmax = Fixed(260),
    lmax = Fixed(90.5),
    linf = Normal(80,10),
    k = Normal(0.6,0.1)
  )),
  
  'Katsuwonus pelamis'
)