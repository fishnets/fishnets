

rm(list=ls())

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
source('collate.R')
source('load_data.R')


# estimate amat in fb

amat_fb <- fb$amat

obj <- function(alpha) lmat - linf * (1 - exp(-k * (alpha - t0)))

for(i in 1:length(amat_fb)) {
  
  cat('row:',i,'\n')
  
 if(is.na(amat_fb[i])) {
   
  lmat <- fb$lmat[i]
  linf <- fb$linf[i]
  k    <- fb$k[i]
  t0   <- fb$t0[i]
  
  if(any(is.na(c(lmat,linf,k,t0)))) next
  if(lmat>linf) next
  if(lmat<(linf * (1 - exp(-k * (0 - t0))))) next
  
  amat_fb[i] <- uniroot(obj,interval=c(0,100))$root
  
 }
  
}



