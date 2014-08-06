
rm(list=ls())

# Source in the package and data
source('collate.R')
source('load_data.R')

# create brt fishnet

brt14 <- Fishnet(
  
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  lmat      = TaxonomicImputer('lmat',c(log,exp)),
  
  linf      = Glmer(log(linf) ~ log(lmax) + log(lmat),exp),
  
  k         = Brter(log(k)~family+trophic+log(linf)+log(amax),exp),
  m         = Brter(log(m)~family+log(linf)+log(k)+log(amax),exp)
  
)
  
names(gs)[match('l',names(gs))]<-'lmat'
gs <- subset(gs,m<2.5)
brt14$fit(gs)

# perform jacknife for different nodes
jknife <- list()
jknife[['linf']] <- brt14$nodes[['linf']]$cross(gs,jacknife=T)
jknife[['k']] <- brt14$nodes[['k']]$cross(gs,jacknife=T)
jknife[['m']]    <- brt14$nodes[['m']]$cross(gs,jacknife=T)

dfr <- lapply(jknife,function(x) x$folds[,c('hat','obs')])
dfr <- rbind(cbind(par='linf',dfr[['linf']]),cbind(par='k',dfr[['k']]),cbind(par='m',dfr[['m']]))

ggplot(dfr) + 
  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~par,scales="free") + 
  theme_bw(base_size=20)


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










