#' A `Node` for matural mortality based on
#' [Hoenig 1983]() p.899
MHoenig1983 <- function(){
  self <- extend(Node,'MHoenig1983')

  self$predictand <- 'm'
  self$predictors <- 'amax'
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,exp(1.46 - 1.01*log(amax)))
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
    
  self$tests <- function(){
  }
  
  self
}
