#' A `Node` for matural mortality based on
#' [Roff 1983]()
MRoff1984 <- function(){
  self <- extend(Node,'MRoff1984')

  self$predictand <- 'm'
  self$predictors <- c('k','lmat','linf')
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,3*k*linf*(1-lmat/linf)/lmat)
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
    
  self$tests <- function(){
  }
  
  self
}
