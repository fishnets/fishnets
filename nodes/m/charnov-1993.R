#' A `Node` for matural mortality based on
#' [Charnov 1993]()
MCharnov1993 <- function(){
  self <- extend(Node,'MCharnov1993')

  self$predictand <- 'm'
  self$predictors <- 'k'
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,1.6*k)
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
    
  self$tests <- function(){
  }
  
  self
}
