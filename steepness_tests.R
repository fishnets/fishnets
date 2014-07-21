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

# geometric mean
gmean <- function(x) exp(mean(log(x),na.rm=T))

# reduce dataset; gometric means for paramters by species
steep_red <- steep_merged %.% 
  select(species, genus, family, class, order, mean_R_z, linf, m, fecundity, trophic, lmat, lmax , k, amax, habit, trophic, depthmax) %.% 
  group_by(order,class,genus,family,species) %.% 
  summarise(mean_R_z = unique(mean_R_z),
            habit = unique(habit),
            trophic = gmean(trophic), 
            linf = gmean(linf), 
            m = gmean(m), 
            depthmax = gmean(depthmax),
            fecundity = gmean(fecundity), 
            trophic = gmean(trophic), 
            lmax = gmean(lmax), 
            lmat = gmean(lmat), 
            k = gmean(k), 
            amax = gmean(amax))


# build a net for steepnes. Use Bayesian nodes in an attempt to not overfit (i.e., to egt better predictive power)

steep_net <- Fishnet(
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
  mean_R_z  = Brter(log(mean_R_z-0.2) ~  habit + log(linf) + log(k) + log(m)+ log(fecundity) +recsigma + trophic+log(depthmax),transform = function(x){exp(x)+0.2},ntrees =15000)
  
)


steep_net$fit(steep_red,impute = T)

# imputed data for all nodes, with only original steepness values
testset <- cbind(subset(steep_net$data,select = -mean_R_z),steep_red['mean_R_z'])

# data with only imputed steepness values, no originals
trainset <- steep_net$data
trainset$mean_R_z[!is.na(steep_red$mean_R_z)] <- NA

#steep_net$nodes$mean_R_z = Bayser(log(mean_R_z-0.2)~log(linf)+log(m)+log(fecundity)+trophic+log(k)+log(depthmax),transform = function(x){exp(x)+0.2})
#steep_net$nodes$mean_R_z$fit(testset)

plot(steep_net$nodes$mean_R_z$brt)
summary(steep_net$nodes$mean_R_z$brt)

Predicted <- steep_net$nodes$mean_R_z$predict(trainset)[!is.na(steep_red$mean_R_z)]
Observed <- steep_red$mean_R_z[!is.na(steep_red$mean_R_z)]

self_pred <- lm(Observed~Predicted)
summary(self_pred)

plot(Predicted,Observed,pch=16)
abline(self_pred$coeff[1],self_pred$coeff[2],col=2,lwd=2)
abline(a=0,b=1,lwd=2)

# not too bad but it's cheating a little
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

steep_cv <- jacknife_cv(testset,steep_net,'mean_R_z')

lm_pred_steep <- lm(Observed~Predicted,data=steep_cv)
summary(lm_pred_steep)

plot(steep_cv,pch=16)
abline(lm_pred_steep$coeff[1],lm_pred_steep$coeff[2],col=2,lwd=2)
abline(0,1,lwd=2)



# check it ------

canary <- steep_net$sample(list(
  species = 'Sebastes pinniger'
),samples = 1000)

ggplot(canary) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,6)) + 
  labs(x='Steepness (z)',y='Density')

salmo <- steep_net$sample(list(
  species = 'Salmo salar'
),samples = 1000)

ggplot(salmo) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,6)) + 
  labs(x='Steepness (z)',y='Density')

Cod <- steep_net$sample(list(
  species = 'Gadus morhua'
),samples = 1000)

ggplot(Cod) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,5)) + 
  labs(x='Steepness (z)',y='Density')


Herring <- steep_net$sample(list(
  species = 'Clupea Harengus'
),samples = 1000)

ggplot(Herring) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,6)) + 
  labs(x='Steepness (z)',y='Density')

####### Bluenose --------

bwa <- steep_net$sample(list(
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
),samples = 10000)

ggplot(bwa) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,4)) + 
  labs(x='Steepness (z)',y='Density')

ggplot(bwa) + 
  geom_point(aes(x=m,y=mean_R_z),alpha=0.4) + 
  scale_x_log10(breaks=seq(0.1,1.1,0.2)) + 
  scale_y_log10(breaks=seq(0.1,1.1,0.2)) + 
  labs(x='Natural mortality rate (M)',y='Steepness')

#### Cross validation-------

sub_data <- subset(mdata,!mdata[[layer]] %in% spp[cv_group==f])


temp_net$fit(sub_data, impute = impute)

cross_val_net <- function(fishnet = NULL, data = NULL, folds = 10, layer = 'species', nodes = 'all', impute = T){
  
  if (nodes == 'all') nodes <- names(fishnet$nodes)
  # cross validate each node
  pred <- vector("list", length(nodes))
  names(pred) <- nodes
  for(name in nodes){
    cat('Cross-validating',name,'\n')    
    # subset data to relevant data
    if (!is.null(fishnet$nodes[[name]]$formula)){
      mdata <- subset(data, id %in% model.frame(paste(paste(deparse(fishnet$nodes[[name]]$formula),collapse=''),'+id'),data)$id)
    } else {mdata <- data}
      
    # random subsets
    spp <- unique(mdata[[layer]])
    cv_group <- sample(1:folds, length(spp), replace = T)
  
    # temp fishnet; don't want to over-write fitting from full net
    temp_net <- fishnet
      
    for (f in 1:folds){
      cat('CV fold',f,'\n')    
      # take out layer data from fold f
      sub_data <- subset(mdata,!mdata[[layer]] %in% spp[cv_group==f])
      temp_net$fit(sub_data, impute = impute)
      rmv <- which(colnames(mdata) == name)
      test_data <- mdata[mdata[[layer]] %in% spp[cv_group==f],-rmv]
      pred[[name]][[f]] <- data.frame(org = mdata[mdata[[layer]] %in% spp[cv_group==f],rmv],pred = temp_net$nodes[[name]]$predict(test_data))
    }
    cat('done','\n') 
  }
}

cross_val_net(fishnet = steep_net, data = steep_net$data, folds = 10, layer = 'species', nodes = 'mean_R_z')

ggplot(cmp) + geom_point(aes(x=preds,y=obs),size=3,alpha=0.3) +
  geom_abline(a=0,b=1) +
  scale_x_log10("Predicted k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2)) + 
  scale_y_log10("Observed k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2))


# Comparison of estiamted vs 'data' -----

steep_net_test <- steep_net

#' Plot density histograms
plot_samples <- function(samples,species_,pars=c('linf','k','m','mean_R_z')){
  data = subset(steep_merged,species==species_)
  
  melted <- melt(samples[,pars])
  data_melted <- melt(data[,pars])
  ggplot(melted,aes(x=value)) +
    geom_bar(data=data_melted,aes(y = ..density..)) +
    geom_density(fill=hsv(0,0.7,0.7),alpha=0.5) +
    facet_wrap(~variable,scales='free') + 
    labs(x='',y='Density') + 
    theme(strip.text.x=element_text(size=10))
}

steep_net_test$fit(subset(steep_merged,species!='Gadus morhua'),impute = T)

preds <- steep_net_test$sample(list(
    species = 'Gadus morhua'
  ))

plot_samples(preds,'Gadus morhua')

plot_samples(
  
  steep_net_test$sample(list(
    species = 'Gadus morhua',
    swimmode = 'subcarangiform',
    habit = 'benthopelagic',
    depthmax = 600,
    temp = 5,
    lmax = 132
  )),
  
  'Gadus morhua'
)

preds <- steep_net_test$sample(dists(
  species =  Fixed('Gadus morhua'),
  swimmode = Fixed('subcarangiform'),
  habit = Fixed('benthopelagic'),
  depthmax = Fixed(600),
  temp = Trapezoid(2,3,6,7),
  lmax = Fixed(132),
  linf = Normal(110,20),
  k = Triangle(0.07,0.13,0.35),
  amax=Fixed(20)
),10000)

plot_samples(preds,'Gadus morhua')

# without imputation

steep_net_noImputation <- steep_net_test
steep_net_noImputation$fit(subset(steep_merged,species!='Gadus morhua'), impute = F)

preds_noImputation <- steep_net_noImputation$sample(dists(
  species =  Fixed('Gadus morhua'),
  swimmode = Fixed('subcarangiform'),
  habit = Fixed('benthopelagic'),
  depthmax = Fixed(600),
  temp = Trapezoid(2,3,6,7),
  lmax = Fixed(132),
  linf = Normal(110,20),
  k = Triangle(0.07,0.13,0.35),
  amax=Fixed(20)
),10000)

plot_samples(preds_noImputation,'Gadus morhua')

###### Skipjack ------

steep_net$fit(subset(steep_merged,species!='Katsuwonus pelamis'))

plot_samples(
  
  steep_net$sample(list(
    species = 'Katsuwonus pelamis',
    family = 'Scombridae'
  )),
  
  'Katsuwonus pelamis'
)

plot_samples(
  
  steep_net$sample(dists(
    species = Fixed('Katsuwonus pelamis'),
    family = Fixed('Scombridae'),
    depthmax = Fixed(260),
    lmax = Fixed(90.5)
  )),
  
  'Katsuwonus pelamis'
)


plot_samples(
  
  steep_net$sample(dists(
    species = Fixed('Katsuwonus pelamis'),
    family = Fixed('Scombridae'),
    depthmax = Fixed(260),
    lmax = Fixed(90.5),
    linf = Normal(80,10),
    k = Normal(0.6,0.1)
  )),
  
  'Katsuwonus pelamis'
)