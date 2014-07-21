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

# Create test net for steepness -----

steep_net <- Fishnet(
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit'),
  depthmax  = TaxonomicImputer('depthmax',c(log,exp)),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp)),
  fecundity = TaxonomicImputer('fecundity',c(log,exp)),
  
  linf      = Glmer(log(linf)~class+order+family+log(lmax),exp),
  k         = Brter(log(k)~class+order+family+log(linf)+habit+log(depthmax)+trophic,exp),
  m         = Svmer(log(m)~class+order+family+log(k)+log(amax),exp),
  lmat      = Glmer(log(lmat)~class+order+family+log(linf),exp),
  #mean_R_z  = Brter(log(mean_R_z-0.2)~log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2},bag.fraction = 0.2,ntrees=0)
  mean_R_z  = Glmer(log(mean_R_z-0.2)~log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2})
  #mean_R_z  = Svmer(log(mean_R_z-0.2)~log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2})
  #recsigma  = RecsigmaThorsonEtAl2014(),
  #recsteep  = RecsteepHeEtAl2006()
  #recauto   = RecautoThorsonEtAl2014()
)
# Fit to Fishbase data -----
steep_net$fit(steep_merged,impute = T)
#steep_net$nodes$mean_R_z$fit(steep_net$data)
#summary(steep_net$nodes$mean_R_z$glm)
summary(steep_net$nodes$linf$glm)
summary(steep_net$nodes$mean_R_z$glm)

testset <- cbind(subset(steep_net$data,select = -mean_R_z),steep_merged['mean_R_z'])
trainset <- steep_net$data

trainset$mean_R_z[!is.na(steep_merged$mean_R_z)] <- NA


steep_net$nodes$mean_R_z = Brter(log(mean_R_z-0.2)~ species+log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2},ntrees=5000)

steep_net$nodes$mean_R_z = Glmer(log(mean_R_z-0.2)~log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2})

steep_net$nodes$mean_R_z$fit(testset)
plot(steep_net$nodes$mean_R_z$glm)
summary(steep_net$nodes$mean_R_z$glm)

# predict back from imputated data to data
steep_net$nodes$mean_R_z$fit(trainset)
pred <- steep_net$nodes$mean_R_z$predict(trainset)[!is.na(steep_merged$mean_R_z)]
org <- steep_merged$mean_R_z[!is.na(steep_merged$mean_R_z)]

plot(pred,org)
abline(lm(org~pred)$coeff[1],lm(org~pred)$coeff[2],col=2,lwd=2)
summary(lm(org~pred))
abline(a=0,b=1,lwd=2)

# how good is imputation?
steep_net$nodes$mean_R_z$fit(steep_merged)
plot(steep_net$nodes$mean_R_z$glm)
summary(steep_net$nodes$mean_R_z$glm)


#summary(steep_net$nodes$mean_R_z$svm)

steep_net_noImpute <- steep_net
steep_net_noImpute$fit(steep_merged,impute = F)
summary(steep_net_noImpute$nodes$linf$glm)

#steep_net$nodes$mean_R_z$fit(steep_merged)
summary(steep_net_noImpute$nodes$mean_R_z$glm)
summary(steep_net_noImpute$nodes$mean_R_z$brt)
summary(steep_net_noImpute$nodes$mean_R_z$svm)

gmean <- function(x) exp(mean(log(x)))

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

steep_net$fit(steep_red,impute = T)

testset <- cbind(subset(steep_net$data,select = -mean_R_z),steep_red['mean_R_z'])
trainset <- steep_net$data

trainset$mean_R_z[!is.na(steep_red$mean_R_z)] <- NA

steep_net$nodes$mean_R_z = Brter(log(mean_R_z-0.2)~ order+log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2},ntrees=5000,bag.fraction=0.8)

steep_net$nodes$mean_R_z = Glmer(log(mean_R_z-0.2)~order+class+log(linf)+log(m)+log(fecundity)+trophic+log(lmat)+log(k)+log(amax),transform = function(x){exp(x)+0.2})

steep_net$nodes$mean_R_z$fit(testset)
plot(steep_net$nodes$mean_R_z$brt,)
summary(steep_net$nodes$mean_R_z$brt)

# predict back from imputated data to data
steep_net$nodes$mean_R_z$fit(trainset)
pred <- steep_net$nodes$mean_R_z$predict(trainset)[!is.na(steep_red$mean_R_z)]
org <- steep_red$mean_R_z[!is.na(steep_red$mean_R_z)]

plot(pred,org)
abline(lm(org~pred)$coeff[1],lm(org~pred)$coeff[2],col=2,lwd=2)
summary(lm(org~pred))
abline(a=0,b=1,lwd=2)

### try a different net

steep_alt <- Fishnet(
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit'),
  depthmax  = TaxonomicImputer('depthmax',c(log,exp)),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp)),
  
  linf      = Bayser(log(linf) ~ f(family,model="iid")+f(class,model="iid")+log(lmax),exp),
  fecundity = Bayser(log(fecundity) ~ f(family,model="iid")+f(class,model="iid") + log(linf) + log(depthmax),exp),
  k         = Bayser(log(k) ~ f(family,model="iid") + log(linf) + f(habit,model="iid") + log(depthmax),exp),
  m         = Bayser(log(m) ~ f(family,model="iid")+f(class,model="iid")+log(k)+log(linf)+f(habit,model="iid")+log(depthmax)+trophic,exp),
  lmat      = Bayser(log(lmat) ~ f(family,model="iid")+log(k)+log(linf)+f(habit,model="iid")+log(depthmax),exp),
  mean_R_z  = Brter(log(mean_R_z-0.2) ~  habit + log(linf) + log(k) + log(m)+ log(fecundity)  + trophic+log(depthmax),transform = function(x){exp(x)+0.2},ntrees =15000)
  
)

steep_alt$fit(steep_red,impute=T)

pred = steep_alt$nodes$m$predict(steep_red)
plot(log(pred),log(steep_red$m))
abline(0,1)

testset <- cbind(subset(steep_alt$data,select = -mean_R_z),steep_red['mean_R_z'])

pred = steep_alt$nodes$mean_R_z$predict(steep_alt$data[!is.na(steep_red$mean_R_z),])
plot((pred),steep_red$mean_R_z[!is.na(steep_red$mean_R_z)],pch = as.numeric(steep_red$order[!is.na(steep_red$mean_R_z)]))
abline(0,1)

# check it ------

canary <- steep_net$sample(list(
  species = 'Sebastes pinniger'
),samples = 10000)

ggplot(canary) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,3)) + 
  labs(x='Steepness (z)',y='Density')

salmo <- steep_net$sample(list(
  species = 'Salmo salar'
),samples = 10000)

ggplot(salmo) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,2)) + 
  labs(x='Steepness (z)',y='Density')

Cod <- steep_net$sample(list(
  species = 'Gadus morhua'
),samples = 10000)

ggplot(Cod) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,5)) + 
  labs(x='Steepness (z)',y='Density')


Herring <- steep_net$sample(list(
  species = 'Clupea Harengus'
),samples = 10000)

ggplot(Herring) + 
  geom_bar(aes(x=mean_R_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,3)) + 
  labs(x='Steepness (z)',y='Density')


ggplot(Cod) + 
  geom_bar(aes(x=m,y=..density..),fill='grey40') + 
  scale_x_continuous() + 
  labs(x='Natural mortality (m)',y='Density')

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