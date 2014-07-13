# steepness tests;

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


steppa <- Steepness$create('./data/steepness')
steppa_merged <- Steepness$merge(steppa,fb)

# Create test net for steepness
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
  mean_R_z = Brter(log(mean_R_z)~class+order+family+log(linf)+m+trophic+lmat+k+amax,transform = exp, bag.fraction = 0.8)
  
  #recsigma  = RecsigmaThorsonEtAl2014(),
  #recauto   = RecautoThorsonEtAl2014()
  #recsteep  = RecsteepHeEtAl2006()
)
# Fit to Fishbase data
steep_net$fit(steppa_merged)

# check it ------


bwa <- steep_net$sample(list(
  species = 'Clupea harengus'
),samples = 1000)

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
