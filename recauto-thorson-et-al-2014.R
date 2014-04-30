#' A `Node` for recruitment frst order autocorrelation (`recauto`) based on
#' [Thorson et al 2014](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0645)
require(truncnorm)

RecautoThorsonEtAl2014 <- function(){
  self <- extend(Node,'RecautoThorsonEtAl2014')
  
  self$predictand = 'recauto'
  self$predictors = 'order'
  
  self$table <- read.csv(text="
order,               mean,  sd
Aulopiformes,       0.464,  0.264
Clupeiformes,       0.435,  0.266
Gadiformes,         0.404,  0.265
Perciformes,        0.466,  0.26
Pleuronectiformes,  0.437,  0.265
Salmoniformes,      0.371,  0.276
Scorpaeniformes,    0.439,  0.264
<other>	            0.426,  0.275",
  stringsAsFactors=F)
  
  self$lookup <- function(data){
    # Look up row in table (in a vectorised way!)
    matches <- match(as.character(data$order),self$table$order)
    matches[is.na(matches)] <- nrow(self$table)
    self$table[matches,]
  }
  
  self$predict <- function(data){
    # Look up mean recsigma in the table
    self$lookup(data)$mean
  }
  
  self$sample <- function(data,samples=1){
    # Lookup the mean and sd in the table
    # Convert these to the scale and shape parameters of a truncated normal distribution
    parameters <- self$lookup(data)
    rtruncnorm(nrow(parameters),a=-0.99,b=0.99,mean=parameters$mean,sd=parameters$sd)
  }
  
  self$tests <- function(){
    # A few simnple tests
    self$predict(list(order='Salmoniformes'))==0.371
    self$predict(list(order='some-unknown-order'))==0.426
    hist(self$sample(list(order='Perciformes'),1000),breaks=100)
  }
  
  self
}
