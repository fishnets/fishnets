#' This is an example using Fishnets for the dark ghost shark (Hydrolagus novaezealandiae)
#' as suggested by Charles Edwards, NIWA. Currently this just examines alternative nodes for 
#' predicting M and does not try to fit a complete Fishnet. It fits a couple of predictor nodes for M
#' from the literature. These will eventually be given there own node class (e.g. a CharnovEtAl2013 node for M which
#' assumes cerain parameter values rather than being fitted)

# Run this from the fishnets main directory
source('../collate.R')

# Load the Fishbase data
fb <- fishbase2000$read('../data/fishbase-2000')

# Fit a GLM representing Charnov et al 2013 equation 3 but with coefficient on `k` allowed
# to differ from 1
# Note that we use `lmat` (length at first maturity) instead of `l` (length of fish at some age in life)
# which is quite fundamentally different from Charnov et al which is as much a model of changes in M
# over life span as it is across species
cea13 <- Glmer(log(m)~log(lmat/linf)+log(k),exp)
cea13$fit(fb)
# Extract coefficients
# They are very different from Charnov's estimates - different data, use of `lmat`?
summary(cea13$glm)

# Fit a GLM representing Gislason et al 2010 equation
# Again note that `lmat` is used instead of `l`. See note above.
gea10 <- Glmer(log(m)~log(lmat)+log(linf)+log(k),exp)
gea10$fit(fb)
# Extract coefficients
# They are different from Gislason's estimates - different data, use of `lmat`?
summary(gea10$glm)

# Fit a more adhoc GLM which allows for different intercepts
# for taxomic groups
glm <- Glmer(log(m)~log(linf)+log(k)+class+order+family,exp)
glm$fit(fb)

# Fit a support vector machine with a bunch of stuff thrown in
svm <- Svmer(log(m)~class+order+family+log(linf)+log(k)+foodtroph,exp)
svm$fit(fb)

# Compile some dark ghost shark data
gsh <- data.frame(
  species='Hydrolagus novaezealandiae',
  genus='Hydrolagus',
  family='Chimaeridae',
  order='Chimaeriformes',
  class='Holocephali'
)
# Estimates from plenary report May 2013
# maturity (cms)
lcrit.male   <- 52.5 
lcrit.female <- 62.5
gsh <- c(gsh,
  lmat = mean(c(lcrit.male,lcrit.female))
)
# Use mean of estimates from Fishbase for other variables
# (since not using an entire Fishnet to populate those)
fbgsh <- subset(fb,family=='Chimaeridae')
#fbgsh <- subset(fb,order=='Chimaeriformes')
#fbgsh <- subset(fb,class=='Holocephali')
gsh <- c(gsh,
   linf = mean(fbgsh$linf),
   k = mean(fbgsh$k),
   foodtroph = mean(fbgsh$foodtroph)
)

# For each of the above predictor nodes for M create a
# distribution of M
nodes = list(cea13=cea13,gea10=gea10,glm=glm,svm=svm)
dists = NULL
for(name in names(nodes)){
  node = nodes[[name]]
  sample = node$sample(gsh,10000)
  dists = rbind(dists,data.frame(node=name,M=sample))
}

# Plot em!
require(ggplot2)
ggplot(dists,aes(x=M,colour=node)) + geom_density() + scale_x_continuous(trans='log')

