#' A `Node` for matural mortality based on
#' [Quinn & Deriso 1999]() 
MQuinn1999 <- function(){
  self <- extend(Node,'MQuinn1999')

  self$predictand <- 'm'
  self$predictors <- 'amax'
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,3/amax)
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
  
  self$tests <- function(){
  }
  
  self
}
