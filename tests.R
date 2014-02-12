#' Some ad-hoc tests and examples...
#' 
#' @author Nokome Bentley

source('collate.R')

# Other packages required
require(ggplot2)
require(reshape)

# Load the fishbase data
fb <- fishbase2000$read('data/fishbase-2000')

# In these tests, examine out-of-sample predictive ability
# for cod
fbcod <- subset(fb,species=='Gadus morhua')
fbnocod <- subset(fb,species!='Gadus morhua')

# Use the various types of Taxon.lookupper nodes
ge <- Genus.lookupper()
ge$fit(fb)
ge$predict(list(
  species='Gadus morhua'
))

fa <- Family.lookupper()
fa$fit(fb)
fa$predict(list(
  genus='Gadus'
))

or <- Order.lookupper()
or$fit(fb)
or$predict(list(
  family='Gadidae'
))

cl <- Class.lookupper()
cl$fit(fb)
cl$predict(list(
  order='Gadiformes'
))

# Test some nodes that predict k....
# Record samples for comparison
ksamples <- list()

ksamples[['fb']] <- sample(fbcod$k,1000,replace=T)

# Create an input list for testing various nodes
cod <- data.frame(
  species='Gadus morhua',
  genus='Gadus',
  family='Gadidae',
  order='Gadiformes',
  class='Actinopterygii'
)

# Taxonomic imputer node
ti <- Taxonomic.imputer('k',c(log,exp),10)
#ti$cross(10,fb)
ti$fit(fbnocod)
ti$predict(cod)
ksamples[['ti']] <- ti$sample(cod,1000)

# The following nodes all predict k but need
# an estimate of linf, so use the mean of observations for cod
cod <- c(cod,
  linf = mean(fbcod$linf,na.rm=T)
)

# Glmer node
glm <- Glmer(log(k)~log(linf)+class+order+family+genus+species,exp,5,5)
#glm$cross(10,fb)
glm$fit(fbnocod)
glm$predict(cod)
ksamples[['glm']] <- glm$sample(cod,1000)

# The following node uses foodtroph, so use the mean of observations for cod
cod <- c(cod,
  foodtroph = mean(fbcod$foodtroph,na.rm=T)
)

# Svmer node
svm <- Svmer(log(k)~log(linf)+class+order+family+genus+species+foodtroph,exp)
#svm$cross(10,fb)
svm$fit(fbnocod)
svm$predict(cod)
ksamples[['svm']] <- svm$sample(cod,1000)

# Compare distributions of k produced by the different types of nodes
ksamples.df <- melt(ldply(ksamples))
names(ksamples.df)[1] <- c('type')
ksamples.df$type = factor(ksamples.df$type,
  levels=c('ti','glm','svm','fb'),
  labels=c('TI','GLM','SVM','G.morhua')
)
ggplot(ksamples.df) + 
  geom_bar(aes(x=value)) +
  scale_x_continuous(trans='log',breaks=c(0.1,0.2,0.3,0.4)) +
  facet_grid(type~.) + 
  labs(x='k',y='Samples')

# Create an example network...

# An example Fishnet using a list of nodes 
fn <- Fishnet(list(
  species = Species.any()
  ,genus = Genus.lookupper()
  ,family = Family.lookupper()
  ,order = Order.lookupper()
  ,class = Class.lookupper()
  
  ,lmax = Taxonomic.imputer('lmax',c(log,exp),3)
  ,tmax = Glmer(log(tmax)~log(lmax)+genus+family+order+class,exp,3,3)
  
  ,linf = Glmer(log(linf)~log(lmax)+species+genus+family+order+class,exp,10,10)
  ,lmat = Glmer(log(lmat)~log(linf)+species+genus+family+order+class,exp,10,10)
  ,k =    Svmer(log(k)~log(linf)+species+genus+family+order+class,exp)
  
  ,m =    Svmer(log(m)~log(tmax)+log(k)+genus+family+order+class,exp)
))

# Fit the entire network without cod data
fn$fit(fbnocod)

# Predict for a generic fish
fn$predict(list(fish=1))

# Predict and sample for a cod
fn$predict(cod)
fn$sample(cod)

# Usually, rather than using fixed variables you want to specify
# prior distributions. So, set up a list of distributions based on
# actual data for cod
cod <- function(which){
  dists <- with(fbcod,list(
    species = Fixed('Gadus morhua'),
    genus = Fixed('Gadus'),
    lmax = Fixed(max(lmax,na.rm=T)),
    tmax = Lognormal(
      log(median(tmax,na.rm=T)),
      0.2
    ),
    linf = Normal(
      mean(linf,na.rm=T),
      sd(linf,na.rm=T)
    ),
    k = Triangle(
      min(k,na.rm=T)*0.9,
      median(k,na.rm=T),
      max(k,na.rm=T)*1.1
    ),
    lmat = Trapezoid(
      min(lmat,na.rm=T)*0.9,
      median(lmat,na.rm=T)*0.9,
      median(lmat,na.rm=T)*1.1,
      max(lmat,na.rm=T)*1.1
    )
  ))[which]
  class(dists) <- c('Distributions','list')
  dists
}

# Generate samples for each level of adding data
samples <- list()
variables <- c('species','genus')
for(add in c('lmax','tmax','linf','lmat','k')){
  variables <- c(variables,add)
  samples[[add]] <- fn$sample(cod(variables),100)
}

# Compare distributions of m produced as additional data is added
samples.df <- ldply(samples)
names(samples.df)[1] <- c('add')
samples.df$add = factor(samples.df$add,
  levels=c("lmax","tmax","linf","lmat","k")
)
ggplot(samples.df) + 
  geom_bar(aes(x=m)) +
  scale_x_continuous(trans='log',breaks=c(0.1,0.2,0.3,0.4)) +
  facet_grid(add~.) + 
  labs(x='Natural mortality (m)',y='Samples')

# Generate a graph of the network
fn$graph("test")

# Store and restore the network
fn$store(list(fish=1),100,'test')
fn$restore('test')
