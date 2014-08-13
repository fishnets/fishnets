#' A `Node` for matural mortality based on
#' [Charnov & Berrigan 1990]()
MCharnov1990 <- function(){
  self <- extend(Node,'MCharnov1990')

  self$predictand <- 'm'
  self$predictors <- 'amat'
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,1.8/amat)
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
    
  self$tests <- function(){
  }
  
  self
}
