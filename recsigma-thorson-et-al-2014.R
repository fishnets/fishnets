#' A `Node` for recruitment variability (`recsigma`) based on
#' [Thorson et al 2014](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0645)

RecsigmaThorsonEtAl2014 <- function(){
  self <- extend(Node,'RecsigmaThorsonEtAl2014')
  
  self$predictand = 'recsigma'
  self$predictors = 'order'
  
  self$table <- read.csv(text="
order,               mean,  sd
Aulopiformes,       0.670, 0.298
Clupeiformes,       0.766, 0.305
Gadiformes,         0.748, 0.293
Perciformes,        0.777, 0.313 
Pleuronectiformes,  0.636, 0.263
Salmoniformes,      0.711, 0.282
Scorpaeniformes,    0.778, 0.318
<other>,            0.737, 0.353",
  stringsAsFactors=F)
  
  self$lookup <- function(data){
    # Look up row in table
    order_ <- as.character(data$order)
    if(order_ %in% self$table$order) subset(self$table,order==order_)
    else subset(self$table,order=='<other>')
  }
  
  self$fit <- function(data){
    # Nothing to be done here
  }
  
  self$predict <- function(data){
    # Look up mean recsigma in the table
    self$lookup(data)$mean
  }
  
  self$sample <- function(data,samples=1){
    # Lookup the mean and sd in the table
    # Convert these to the scale and shape parameters of a Gamma distribution
    with(self$lookup(data),{
      shape <- (mean/sd)^2
      scale <- (sd^2)/mean
      rgamma(samples,scale=scale,shape=shape)
    })
  }
  
  self$tests <- function(){
    # A few simnple tests
    self$predict(list(order='Gadiformes'))==0.748
    self$predict(list(order='some-unknown-order'))==0.737
    hist(self$sample(list(order='Perciformes'),1000),breaks=100)
  }
  
  self
}
