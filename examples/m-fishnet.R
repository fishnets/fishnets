
rm(list=ls())

# Source in the package and data
source('collate.R')
source('load_data.R')

# function to make testset for cross validation (PN)
make_testset <- function(net,data,name){
  testset <- data.frame(net$data[,-which(colnames(net$data) == name)], name = data[name])
  testset
}



# create brt fishnet

brt14 <- Fishnet(
  
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp)),
  
  linf      = Brter(log(linf) ~ class+order+family+log(lmax),exp),
  lmat      = Brter(log(lmat) ~ class+order+family+log(linf),exp),

  m         = Brter(log(m)~family+log(linf)+log(lmat)+log(k)+log(amax),exp)
  
)
  
names(gs)[match('l',names(gs))]<-'lmat'
brt14$fit(gs)

# perform jacknife for different nodes
jknife <- list()
jknife[['linf']] <- brt14$nodes[['linf']]$cross(gs,jacknife=T)
jknife[['lmat']] <- brt14$nodes[['lmat']]$cross(gs,jacknife=T)
jknife[['m']]    <- brt14$nodes[['m']]$cross(gs,jacknife=T)

dfr <- lapply(jknife,function(x) x$folds[,c('hat','obs')])
dfr <- rbind(cbind(par='linf',dfr[['linf']]),cbind(par='lmat',dfr[['lmat']]),cbind(par='m',dfr[['m']]))

ggplot(dfr) + 
  geom_point(aes(x=obs,y=hat),size=3,alpha=0.3) + 
  geom_abline(a=0,b=1) + 
  facet_wrap(~par,scales="free")



# predictions
bwa <- brt14$sample(list(
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
  geom_bar(aes(x=mean_BH_z,y=..density..),fill='grey40') + 
  scale_x_continuous(limits=c(0.2,1)) + 
  labs(x='Steepness (z)',y='Density')

ggplot(bwa) + 
  geom_point(aes(x=m,y=mean_BH_z),alpha=0.4) + 
  scale_x_log10(breaks=seq(0.1,1.1,0.2)) + 
  scale_y_log10(breaks=seq(0.1,1.1,0.2)) + 
  labs(x='Natural mortality',y='Steepness')

