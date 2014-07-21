# steepness tests;

# Source in the package ----
source('collate.R')

# Load the Fishbase data

fb <- FishbaseWeb$read('data/fishbase-web')
# An an id column for indexing later
fb$id <- 1:nrow(fb)
# Add a dummy row for helping with predictor nodes
# that need to have at least two predictors
fb$dummy <- 1.0

# delete some k and linf values for testing purposes

fb$linf[sample(x = nrow(fb),size = 2000)] <- NA
fb$k[sample(x = nrow(fb),size = 2000)] <- NA

# Create test net for imputation absed on Bea14-----

impute_net <- Fishnet(
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
  
  linf      = Bayser(log(linf)~f(family,model='iid')+f(order,model='iid')+log(lmax),exp),
  k         = Bayser(log(k) ~ f(family,model='iid')+log(linf)+habit+log(depthmax)+trophic,exp),
  m         = Bayser(log(m) ~ f(family,model='iid')+log(k)+log(amax),exp),
  lmat      = Glmer(log(lmat)~class+order+family+log(linf),exp)
)
# Fit to Fishbase data -----

no_impute_net <- impute_net
# with imputation
impute_net$fit(fb,impute = T)
# without imputation
no_impute_net$fit(fb)

# Cross validation; doesn't seem to matter much?

cv_impute <- impute_net$nodes$m$cross(fb,20)
cv_impute2 <- impute_net$nodes$m$cross(impute_net$nodes$m$fit_data,20)  

cv_no_impute <- no_impute_net$nodes$m$cross(fb,20)


