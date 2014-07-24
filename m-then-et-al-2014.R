#' A `Node` for matural mortality based on
#' [Then et al 2014]() 
MThenEtAl2014 <- function(){
  self <- extend(Node,'MThenEtAl2014')

  self$predictand <- 'm'
  self$predictors <- c('amax')
  
  self$fit <- function(data){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    rep(1,nrow(data))
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
  
  self$tests <- function(){
  }
  
  self
}
