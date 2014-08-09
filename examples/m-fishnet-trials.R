
rm(list=ls())

# Source in the package and data
source('collate.R')
source('load_data.R')

##############
# groom data #
##############

fb[which(fb$temp<=0),'temp'] <- NA
fb[which(fb$m>2),'m']        <- NA

gs[which(gs$m>2),'m']        <- NA


#############
# fishnet 1 #
#############

mfishnet1 <- Fishnet(
  
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit',c(log,exp)),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  temp      = TaxonomicImputer('temp',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  
  amax      = Brter(log(amax)~family+habit+trophic+log(temp)+log(lmax),exp,ntrees=5000,learning.rate=0.001),
  k         = Brter(log(k)~family+log(temp)+log(lmax)+log(amax),exp,ntrees=5000,learning.rate=0.001),
  m         = Brter(log(m)~family+log(temp)+log(k)+log(amax),exp,ntrees=5000,learning.rate=0.001)
  
)

# perform jacknife for different nodes
jknife <- list()
jknife[['fb']] <- list()
jknife[['fb']][['amax']] <- mfishnet1$nodes[['amax']]$cross(fb,folds=10)
jknife[['fb']][['k']]    <- mfishnet1$nodes[['k']]$cross(fb,folds=10)
jknife[['fb']][['m']]    <- mfishnet1$nodes[['m']]$cross(fb,folds=10)
jknife[['gs']] <- list()
jknife[['gs']][['amax']] <- mfishnet1$nodes[['amax']]$cross(gs,jacknife=T)
jknife[['gs']][['k']]    <- mfishnet1$nodes[['k']]$cross(gs,jacknife=T)
jknife[['gs']][['m']]    <- mfishnet1$nodes[['m']]$cross(gs,jacknife=T)

saver(jknife,name='mfishnet1_jknife')


dfr <- lapply(jknife,function(x) x$folds[,c('hat','obs')])
dfr <- rbind(cbind(par='amax',dfr[['amax']]),cbind(par='k',dfr[['k']]),cbind(par='m',dfr[['m']]))

ggplot(dfr) + 
  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~par,scales="free") + 
  theme_bw(base_size=20)

# cross validate fishnet

data.fb <- fb[,c('m','species','genus','family','order','class','habit','trophic','temp','lmax','k','amax')]
data.fb <- data.fb[!is.na(data.fb$m),]

cvfishnet <- function(data,folds=10,byspecies=F) {
  
  if(byspecies) {
    data$species <- factor(data$species)
    spp <- levels(data$species)
    folds <- length(spp)
    levels(data$species) <- 1:folds
    folder <- as.numeric(data$species)
    levels(data$species) <- spp
  } else {
    folder = rep(1:folds,length.out=nrow(data))
    folder = sample(folder)
  }
  
  # Result data.frame
  results = NULL
  # For each fold...
  for(fold in 1:folds){
    cat("Fold",fold,": ")
    # Define training a testing datasets
    train = data[folder!=fold,]
    test = data[folder==fold,]
    
    # test vector
    tests <- test[,'m']
    test[,'m'] <- NA
    
    # Fit model to training data
    mfishnet1$fit(train)
    
    # Get nodes to predict values for their predictands
    # node prediction 'fills in the blanks' at each step
    # i.e. existent data values are not over-written
    test[,'amax'] <- mfishnet1$nodes[['amax']]$predict.safe(test,na.strict=T,na.keep=T)
    test[,'k']    <- mfishnet1$nodes[['k']]$predict.safe(test,na.strict=T,na.keep=T)
    test[,'m']    <- mfishnet1$nodes[['m']]$predict(test,na.strict=T,na.keep=T)
    
    # get predictions
    preds <- test[,'m']
    
    # remove NA's
    na.loc <- is.na(preds) | is.na(tests)
    preds <- preds[!na.loc]
    tests <- tests[!na.loc]
    
    # Calculate various prediction errors
    if(!all(na.loc)) {
      me = mean(abs(tests-preds))
      mse = mean((tests-preds)^2)
      mpe = mean(abs((tests-preds)/tests))
      r2 = cor(tests,preds,use="pairwise.complete.obs")^2
      dev = mean((tests - preds) * (tests - preds))
      
      # Add to results
      results = rbind(results,data.frame(
        fold = fold,
        obs = mean(tests),
        hat = mean(preds),
        me = me,
        mse = mse,
        mpe = mpe,
        r2 = r2,
        dev = dev
      ))
    }
    
  }
  # Summarise results
  summary <- data.frame(mean=apply(results[,4:8],2,mean,na.rm=T),se=apply(results[,4:8],2,function(x) sd(x)/sqrt(length(x))))
  # Return summary and raw results
  return (list(
    summary = summary,
    folds = results
  ))
}

#mfishnet1.res           <- cvfishnet(data.fb)
mfishnet1.res.byspecies <- cvfishnet(data.fb,byspecies=T)

saver(mfishnet1.res.byspecies,name='mfishnet1.res')

#dfr <- mfishnet1.res$folds[,c('hat','obs')]
dfr <- data.frame(mfishnet1.res.byspecies$folds[,c('hat','obs')],label='Predicted Values')
dfr <- rbind(dfr,data.frame(obs=dfr$obs,hat=dfr$hat-dfr$obs,label='Prediction Residuals'))

fig <- ggplot(dfr,aes(x=obs,y=hat)) + 
  geom_point(size=4,alpha=0.3) + 
  geom_abline(aes(intercept=a,slope=b),data=data.frame(label=c('Predicted Values','Prediction Residuals'),a=c(0,0),b=c(1,0)),colour="#990000", linetype="dashed",size=1) +
  facet_wrap(~label,scale='free_y',ncol=1) +
  theme_bw(base_size=20) +
  labs(x='Observed value',y='')

pdfr(fig,width=12,name='mfishnet1_res_byspecies')


# SOME EXAMPLES
########################################

brt14$fit(fb)

bwa <- brt14$sample(list(
  species = 'Hyperoglyphe antarctica'
),samples = 10000)

ggplot(bwa) + 
  geom_bar(aes(x=m,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.2)) + 
  labs(x='Natural mortality rate (M)',y='Density')

########################################

bwa <- brt14$sample(list(
  species = 'Hyperoglyphe antarctica',
  # Maximum length, temperature and 
  # maximum depth from Fishbase
  lmax = 140,
  #temp = 11,
  #depthmax = 1500,
  # Female growth and max age from 
  # Horn et al 2010
  #linf = 92.5,
  #k = 0.071,
  amax = 71  
),samples = 10000)

ggplot(bwa) + 
  geom_bar(aes(x=m,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.2)) + 
  labs(x='Natural mortality rate (M)',y='Density')

########################################

# OUT OF SAMPLE PREDICTION
###########################################

#' Plot density histograms
plot_samples <- function(samples,species_,pars=c('linf','k','m')){
  data = subset(fb,species==species_)
  
  melted <- melt(samples[,pars])
  data_melted <- melt(data[,pars])
  ggplot(melted,aes(x=value)) +
    geom_bar(data=data_melted,aes(y = ..density..)) +
    geom_density(fill=hsv(0,0.7,0.7),alpha=0.5) +
    facet_wrap(~variable,scales='free') + 
    labs(x='',y='Density') + 
    theme_bw(base_size=20)
}

brt14$fit(subset(fb,species!='Gadus morhua'))

plot_samples(
  
  brt14$sample(list(
    species = 'Gadus morhua',
    lmax = 132
  )),
  
  'Gadus morhua'
)

brt14$fit(subset(fb,species!='Katsuwonus pelamis'))

plot_samples(
  
  brt14$sample(list(
    species = 'Katsuwonus pelamis',
    family = 'Scombridae',
    lmax = 90.5
  )),
  
  'Katsuwonus pelamis'
)










