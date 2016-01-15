#' A `Node` for matural mortality based on
#' [Jensen 2001]()
MJensen2001 <- function(){
  self <- extend(Node,'MJensen2001')

  self$predictand <- 'm'
  self$predictors <- c('k','temp')
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,exp(0.66*log(k)+0.45*log(temp)))
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
    
  self$tests <- function(){
  }
  
  self
}
