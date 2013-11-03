load("data/fishbase-2000/fishbase-2000.RData")
names(fb) = tolower(names(fb))
names(fb)[names(fb)=='lm'] = 'lmat'
names(fb)[names(fb)=='loo'] = 'linf'
fb$species = factor(with(fb,paste(genus,species)))

cl <- Class.lookupper()
cl$fit(fb)
cl$predict(list(order='Clupeiformes'))

ti <- Taxonomic.imputer('lmax',c(log,exp),10)
#ti$cross(5,fb)
ti$fit(fb)
data <- data.frame(
  species='Pagrus auratus',
  genus='Pagrus',
  family='Sparidae',
  order='Perciformes',
  class='Actinopterygii'
)
ti$predict(data)
ti$sample(data,10)
hist(ti$sample(data,10000),breaks=30)

glm <- Glmer(log(k)~log(linf)+species+genus+family+order+class,exp,10,10)
#glm$cross(10,fb)
glm$fit(fb)
data <- data.frame(
  species='Pagrus auratus',
  genus='Pagrus',
  family='Sparidae',
  order='Perciformes',
  class='Actinopterygii',
  linf = 60
)
glm$predict(data)
glm$sample(data,10)
hist(glm$sample(data,1000),breaks=30)

fn <- Fishnet(list(
  genus = Genus.lookupper()
  ,family = Family.lookupper()
  ,order = Order.lookupper()
  ,class = Class.lookupper()
  
  ,lmax = Taxonomic.imputer('lmax',c(log,exp),3)
  ,linf = Glmer(log(linf)~log(lmax)+species+genus+family+order+class,exp,10,10)
  ,lmat = Glmer(log(lmat)~log(linf)+species+genus+family+order+class,exp,10,10)
  ,k =    Glmer(log(k)~log(linf)+species+genus+family+order+class,exp,10,10)
  ,tmax = Glmer(log(tmax)~log(lmax)+genus+family+order+class,exp,3,3)
  ,m =    Glmer(log(m)~log(tmax)+log(k)+genus+family+order+class,exp,3,3)
))

fn$fit(fb)

fn$predict(list(
  species = 'Pagrus auratus',
  genus = 'Pagrus',
  linf = 60,
  k = 0.1,
  tmax = 60
))

temp <-fn$sample(list(
  species = 'Katsuwonus pelamis'
),100)
hist(temp$m,breaks=30)

temp <- fn$sample(dists(
  species = Fixed('Katsuwonus pelamis'),
  linf = Normal(90,10),
  k = Normal(0.4,0.1)
),100)
hist(temp$m,breaks=30)

