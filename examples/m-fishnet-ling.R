
rm(list=ls())

# Source in the package and data
source('collate.R')
source('load_data.R')

# create brt fishnet for ling
# available data
subset(fb,family=='Ophidiidae')

# fishnet
brt14 <- Fishnet(
  
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  temp      = TaxonomicImputer('temp',c(log,exp)),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  sex       = TaxonomicImputer('sex',c(log,exp)),
  
  linf      = Glmer(log(linf) ~ sex + log(lmax),exp),
  
  k         = Brter(log(k)~family+trophic+sex+log(temp)+log(linf)+log(amax),exp),
  m         = Brter(log(m)~family+sex+log(linf)+log(k)+log(amax),exp)
  
)

# PREDICTIONS

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


brt14$fit(subset(fb,species!='Genypterus blacodes'))

brt14$predict(list(
  species = 'Genypterus blacodes',
  family='Ophidiidae',
  amax=30,sex='M'
))

plot_samples(
  
  brt14$sample(list(
    species = 'Genypterus blacodes',
    family='Ophidiidae',
    amax=30,sex='M',
    linf=113.9,k=0.127
  )),
  
  'Genypterus blacodes'
)

plot_samples(
  
  brt14$sample(list(
    species = 'Genypterus blacodes',
    family='Ophidiidae',
    amax=30,sex='F',
    linf=156.4,k=0.083
  )),
  
  'Genypterus blacodes'
)





