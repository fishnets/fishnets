
#' A `Node` for recruitment variability (`recsigma`) based on
#' [Thorson et al 2014](http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2013-0645)
RecsigmaThorsonEtAl2014 <- function(){
  self <- extend(Node,'RecsigmaThorsonEtAl2014')
  
  self$predictors = 'order'
  self$predictand = 'recsigma'
  
  # Table 2 of Thorson et al 2014
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
    # Look up row in table (in a vectorised way!)
    matches <- match(as.character(data$order),self$table$order)
    matches[is.na(matches)] <- nrow(self$table)
    self$table[matches,]
  }
  
  self$predict <- function(data){
    # Look up mean recsigma in the table
    self$lookup(data)$mean
  }
  
  self$sample <- function(data){
    # Lookup the mean and sd in the table and 
    # then convert to the scale and shape parameters
    # of a Gamma distribution
    parameters <- self$lookup(data)
    shape <- with(parameters,(mean/sd)^2)
    scale <- with(parameters,(sd^2)/mean)
    rgamma(nrow(parameters),scale=scale,shape=shape)
  }
  
  self$tests <- function(){
    # A few simple tests
    self$predict(list(order='Gadiformes'))==0.748
    self$predict(list(order='some-unknown-order'))==0.737
    hist(self$sample(list(order='Perciformes'),1000),breaks=100)
  }
  
  self
}
