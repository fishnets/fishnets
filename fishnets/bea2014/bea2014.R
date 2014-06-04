#' Development of a Fishnet for Bentley et al 2014 (hopefully!)
#' 
#' @author Nokome Bentley
require(ggplot2)

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
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

# Create Bea14
Bea14 <- Fishnet(
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
  
  linf      = Glmer(log(linf)~class+order+family+log(lmax),exp),
  k         = Brter(log(k)~class+order+family+log(linf)+habit+log(depthmax)+trophic,exp),
  m         = Svmer(log(m)~class+order+family+log(k)+log(amax),exp),
  lmat      = Glmer(log(lmat)~class+order+family+log(linf),exp),
  
  recsigma  = RecsigmaThorsonEtAl2014(),
  recauto   = RecautoThorsonEtAl2014()
  #recsteep  = RecsteepHeEtAl2006()
)
# Fit to Fishbase data
Bea14$fit(fb)

########################################

bwa <- Bea14$sample(list(
  species = 'Hyperoglyphe antarctica'
),samples = 10000)

ggplot(bwa) + 
  geom_bar(aes(x=m,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.2)) + 
  labs(x='Natural mortality rate (M)',y='Density')

ggplot(bwa) + 
  geom_point(aes(x=m,y=k),alpha=0.4) + 
  scale_x_log10(breaks=seq(0.1,1.1,0.2)) + 
  scale_y_log10(breaks=seq(0.1,1.1,0.2)) + 
  labs(x='Natural mortality rate (M)',y='Growth rate (k)')

########################################

bwa <- Bea14$sample(list(
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
  geom_bar(aes(x=m,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.2)) + 
  labs(x='Natural mortality rate (M)',y='Density')

ggplot(bwa) + 
  geom_point(aes(x=m,y=k),alpha=0.4) + 
  scale_x_log10(breaks=seq(0.1,1.1,0.2),limits=c(0.01,1.5)) + 
  scale_y_log10(breaks=seq(0.1,1.1,0.2),limits=c(0.01,1.5)) + 
  labs(x='Natural mortality rate (M)',y='Growth rate (k)')

########################################

bwa <- Bea14$sample(dists(
  species = Fixed('Hyperoglyphe antarctica'),
  # Maximum length, temperature and 
  # maximum depth from Fishbase
  lmax = Fixed(140),
  temp = Uniform(9,13),
  depthmax = Trapezoid(1200,1300,1400,1500),
  # Female growth and max age from 
  # Horn et al 2010
  linf = Normal(92.5,10),
  k = Triangle(0.05,0.07,0.09),
  amax = Uniform(65,75)
),samples = 10000)

ggplot(bwa) + 
  geom_bar(aes(x=k,y=..density..),fill='grey40') + 
  labs(x='Growth coefficient (k)',y='Density')

ggplot(bwa) + 
  geom_point(aes(x=m,y=k),alpha=0.4) + 
  scale_x_log10(breaks=seq(0.1,1.1,0.2),limits=c(0.01,1.5)) + 
  scale_y_log10(breaks=seq(0.1,1.1,0.2),limits=c(0.01,1.5)) + 
  labs(x='Natural mortality rate (M)',y='Growth rate (k)')

########################################

Fishnet(
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  linf      = Glmer(log(linf)~log(lmax),exp),
  k         = Glmer(log(k)~log(linf),exp),
  m         = Glmer(log(m)~log(k),exp)
)$graph()

########################################

Bea14$graph()

########################################

# Example of TaxonomicImputer
rs <- Fishnet(
  k = TaxonomicImputer('k',c(log,exp),10)
)$fit(fb)$sample(list(
  species = 
),samples = 10000)

ggplot(rs) + 
  geom_density(aes(x=recsigma,colour=order)) + 
  labs(x='Stock-recruitment variability (recsigma)',y='Density',colour="") + 
  theme(legend.position='top')

########################################

fb_full = subset(fb,id %in% model.frame(
  log(k) ~ id + class + order + family + genus + species + 
    swimmode + habit + feeding + diet + trophic + log(depthmax) + temp + 
    fecundity + log(lmax)
,fb)$id)
nrow(fb_full)
length(unique(fb_full$species))
fb_full$species = factor(fb_full$species)
spp <- unique(fb_full$species)
spp_group <- sample(1:2,length(spp),replace=T)
sp1 <- subset(fb_full,species %in% spp[spp_group==1])
sp2 <- subset(fb_full,species %in% spp[spp_group==2])

brt <- 

Brter(
  log(k) ~ class + order + family + genus + species + 
     swimmode + habit + feeding + diet + trophic + log(depthmax) + temp + 
     fecundity + log(lmax),
exp,ntrees=2000)

brt$fit(fb_full)
cmp <- data.frame(
  preds = brt$predict(fb_full),
  obs = fb_full$k
)
ggplot(cmp) + geom_point(aes(x=preds,y=obs),size=3,alpha=0.3) +
  geom_abline(a=0,b=1) +
  scale_x_log10("Predicted k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2)) + 
  scale_y_log10("Observed k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2))


par(las=1,mar=c(4,7,1,1))
summary(brt$brt)

brt$fit(sp1)
cmp <- data.frame(
  preds = brt$predict(sp2),
  obs = sp2$k
)
ggplot(cmp) + geom_point(aes(x=preds,y=obs),size=3,alpha=0.3) +
  geom_abline(a=0,b=1) +
  scale_x_log10("Predicted k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2)) + 
  scale_y_log10("Observed k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2))

brt <- 

Brter(
  log(k) ~ class + order + 
    swimmode + habit + feeding + diet + trophic + log(depthmax) + temp + 
    fecundity + log(lmax),
  exp,ntrees=2000)


brt$fit(fb_full)
par(las=1,mar=c(4,7,1,1))
summary(brt$brt)

brt$fit(sp1)
cmp <- data.frame(
  preds = brt$predict(sp2),
  obs = sp2$k
)
ggplot(cmp) + geom_point(aes(x=preds,y=obs),size=3,alpha=0.3) +
  geom_abline(a=0,b=1) +
  scale_x_log10("Predicted k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2)) + 
  scale_y_log10("Observed k",breaks=c(0.1,0.2,0.5,1.0,2.0),limits=c(0.05,2))

########################################

# Get previously generated results to plot
load('fishnets/bea2014/bea2014.RData')
ggplot(summary) + 
  geom_point(aes(x=term,y=r2,colour=method,shape=method),size=3) + 
  ylim(0,1) + 
  scale_shape_manual(values=1:10) + 
  labs(x="Term",y="Cross validation correlation",colour="Method",shape="Method") + 
  theme(axis.text=element_text(angle=90))


########################################

# Example of recsigma priors from Thorson et al
rs <- Fishnet(
  recsigma = RecsigmaThorsonEtAl2014()
)$sample(list(
  order = c('Salmoniformes','Pleuronectiformes','Scorpaeniformes')
),samples = 10000)

ggplot(rs) + 
  geom_density(aes(x=recsigma,colour=order)) + 
  labs(x='Stock-recruitment variability (recsigma)',y='Density',colour="") + 
  theme(legend.position='top')

###########################################
# Comparison of 

#' Plot density histograms
plot_samples <- function(samples,species_,pars=c('linf','k','m','lmat')){
  data = subset(fb,species==species_)
  
  melted <- melt(samples[,pars])
  data_melted <- melt(data[,pars])
  ggplot(melted,aes(x=value)) +
    geom_bar(data=data_melted,aes(y = ..density..)) +
    geom_density(fill=hsv(0,0.7,0.7),alpha=0.5) +
    facet_wrap(~variable,scales='free') + 
    labs(x='',y='Density') + 
    theme(strip.text.x=element_text(size=10))
}

Bea14$fit(subset(fb,species!='Gadus morhua'))

plot_samples(
  
  Bea14$sample(list(
    species = 'Gadus morhua'
  )),
  
  'Gadus morhua'
)

plot_samples(
  
  Bea14$sample(list(
    species = 'Gadus morhua',
    swimmode = 'subcarangiform',
    habit = 'benthopelagic',
    depthmax = 600,
    temp = 5,
    lmax = 132
  )),
  
  'Gadus morhua'
)

plot_samples(
  
  Bea14$sample(dists(
    species =  Fixed('Gadus morhua'),
    swimmode = Fixed('subcarangiform'),
    habit = Fixed('benthopelagic'),
    depthmax = Fixed(600),
    temp = Trapezoid(2,3,6,7),
    lmax = Fixed(132),
    linf = Normal(110,20),
    k = Triangle(0.07,0.13,0.35)
  )),
  
  'Gadus morhua'
)

###############################

Bea14$fit(subset(fb,species!='Katsuwonus pelamis'))

plot_samples(
  
  Bea14$sample(list(
    species = 'Katsuwonus pelamis',
    family = 'Scombridae'
  )),
  
  'Katsuwonus pelamis'
)

plot_samples(
  
  Bea14$sample(dists(
    species = Fixed('Katsuwonus pelamis'),
    family = Fixed('Scombridae'),
    depthmax = Fixed(260),
    lmax = Fixed(90.5)
  )),
  
  'Katsuwonus pelamis'
)


plot_samples(
  
  Bea14$sample(dists(
    species = Fixed('Katsuwonus pelamis'),
    family = Fixed('Scombridae'),
    depthmax = Fixed(260),
    lmax = Fixed(90.5),
    linf = Normal(80,10),
    k = Normal(0.6,0.1)
  )),
  
  'Katsuwonus pelamis',
  pars = c('linf','k','m','lmat')
)
