
rm(list=ls())

require(ggplot2)

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
source('collate.R')

# Load the Fishbase data
fb <- FishbaseWeb$read('data/fishbase-web')
# Load Gislason data
gs <- GislasonEtAl2010Data$read('data/gislason-et-al-2010')

# Create a fishnet that can be used to impute values
imputer <- Fishnet(
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit'),
  depthmax  = TaxonomicImputer('depthmax',c(log,exp)),
  temp      = TaxonomicImputer('temp'),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp))
)
imputer$fit(fb)
gs <- imputer$predict(gs)


# Source in m related nodes (they might nt be in collate.R yet)
# charnov 2013
source('m-charnov-et-al-2013.R')
source('m-charnov-et-al-2013-fitted.R')

cea13 <- MCharnovEtAl2013()

with(fb,plot(log(m/k)~log(lmat/linf)))
cea13fit <- MCharnovEtAl2013Fitted()
cea13fit$fit(fb)

cea13$cross(fb)
cea13fit$cross(fb)

with(gs,plot(log(m/k)~log(l/linf)))
gs$lmat = gs$l
lm(log(m/k)~log(l/linf),data=gs)

cea13fit$cross(fb)

# hoenig 1983
source('m-hoenig-1983.R')
source('m-hoenig-1983-fitted.R')

h83 <- MHoenig1983()
h83fit <- MHoenig1983Fitted()

hist(h83$predict(fb))
h83fit$fit(fb)
hist(h83fit$predict(fb))

h83$cross(fb)

h83fit$cross(fb)

# boosted regression trees
brt.min <- Brter(log(m)~log(k)+log(amax),exp)
fb.min <- subset(fb,!is.na(k) & !is.na(amax))
brt.min$fit(fb.min)
brt.min$cross(fb.min)
hist(brt.min$sample(brt.min$expand(data.frame(k=.3,amax=5),1000)))
summary(brt.min$brt)
nrow(fb.min)

brt.full <- Brter(log(m)~class+order+family+log(linf)+log(temp)+log(k)+log(amax),exp)
fb.full <- subset(fb,!is.na(class) & !is.na(order) & !is.na(family) & !is.na(linf) & !is.na(temp) & temp>0 & !is.na(k) & !is.na(amax))
brt.full$fit(fb.full)
brt.full$cross(fb.full)
hist(brt.full$sample(brt.full$expand(data.frame(k=.3,amax=5),1000)))
summary(brt.full$brt)
nrow(fb.full)

