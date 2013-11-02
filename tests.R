load("data/fishbase-2000/fishbase-2000.RData")
names(fb) = tolower(names(fb))
names(fb)[names(fb)=='lm'] = 'lmat'
names(fb)[names(fb)=='loo'] = 'linf'
fb$species = with(fb,paste(genus,species))

cl <- Class.lookupper()
cl$fit(fb)
cl$predict(list(order='Clupeiformes'))

ti <- Taxonomic.imputer('k',3)
ti$cross(5,fb)
ti$fit(fb)
ti$predict(data.frame(
  species='Pagrus auratus',
  genus='Pagrus',
  family='Sparidae',
  order='Perciformes',
  class='Actinopterygii'
))

glm <- Glmer(log(k)~log(linf)+order+family,exp,3,3)
glm$cross(10,fb)
glm$fit(fb)
glm$predict(list(
  order = 'Perciformes',
  linf = 40
))

fn <- Fishnet(list(
  genus = Genus.lookupper()
  ,family = Family.lookupper()
  ,order = Order.lookupper()
  ,class = Class.lookupper()
  
  ,lmax = Taxonomic.imputer('lmax',3)
  ,linf = Glmer(log(linf)~log(lmax),exp)
  ,lmat = Glmer(log(lmat)~log(linf)+order,exp)
  ,k = Glmer(log(k)~log(linf)+order,exp)
  ,tmax = Glmer(log(tmax)~log(lmax)+order,exp)
  ,m = Glmer(log(m)~log(tmax)+log(k)+order,exp)
))

fn$fit(fb)

fn$predict(
  species = 'Pagrus auratus',
  genus = 'Pagrus',
  linf = 60,
  k = 0.1,
  tmax = 60
)

fn$predict(
  species = 'Katsuwonus pelamis'
)

