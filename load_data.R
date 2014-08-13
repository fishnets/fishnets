

# Load the Fishbase data
fb <- FishbaseWeb$read('data/fishbase-web')
# Load Gislason data
gs <- GislasonEtAl2010Data$read('data/gislason-et-al-2010')

# Create a fishnet that can be used to impute columns that 
# are missing in the Gislason data. This fishnet is intended to be 
# simple and mostly uses TaxonomicImputer
imputer <- Fishnet(
  species   = SpeciesRandom(),
  genus     = GenusParser(),
  family    = FamilyLookupper(),
  order     = OrderLookupper(),
  class     = ClassLookupper(),
  
  habit     = TaxonomicImputer('habit'),
  depthmax  = TaxonomicImputer('depthmax',c(log,exp)),
  trophic   = TaxonomicImputer('trophic',c(log,exp)),
  lmat      = TaxonomicImputer('lmat',c(log,exp)),
  lmax      = TaxonomicImputer('lmax',c(log,exp)),
  amax      = TaxonomicImputer('amax',c(log,exp))
)
# Do imputation base on Fishbase data
imputer$fit(fb)
gs <- imputer$predict(gs)
