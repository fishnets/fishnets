
##############
# load data  #
##############

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
# Do imputation based on Fishbase data
imputer$fit(fb)
gs <- imputer$predict(gs)
rm(imputer)


##############
# groom data #
##############

fb[which(fb$temp<=0),'temp'] <- NA
fb[which(fb$m>1.5),'m']        <- NA

gs[which(gs$m>1.5),'m']        <- NA

# calculate amat in fishbase using VB growth equation
obj <- function(alpha) lmat - linf * (1 - exp(-k * (alpha - t0)))

for(i in 1:length(fb$amat)) {
  
  if(is.na(fb$amat[i])) {
    
    lmat <- fb$lmat[i]
    linf <- fb$linf[i]
    k    <- fb$k[i]
    t0   <- fb$t0[i]
    
    if(any(is.na(c(lmat,linf,k,t0)))) next
    if(lmat>linf) next
    if(lmat<(linf * (1 - exp(-k * (0 - t0))))) next
    
    fb$amat[i] <- uniroot(obj,interval=c(0,100))$root
    
  }
}

# estimate amat in gislasson database using VB growth
# equation and assuming t0 = 0 (Gislason 2010)
gs$amat <- NA

for(i in 1:length(gs$amat)) {
  
  if(is.na(gs$amat[i])) {
    
    lmat <- gs$lmat[i]
    linf <- gs$linf[i]
    k    <- gs$k[i]
    t0   <- 0
    
    if(any(is.na(c(lmat,linf,k,t0)))) next
    if(lmat>linf) next
    if(lmat<(linf * (1 - exp(-k * (0 - t0))))) next
    
    gs$amat[i] <- uniroot(obj,interval=c(0,100))$root
    
  }
}

rm(lmat,linf,k,t0,i,obj)

save(fb,file='fb_data.Rdata')
save(gs,file='gs_data.Rdata')
