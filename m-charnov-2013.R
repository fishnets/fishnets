#' A `Node` for matural mortality based on
#' [Charnov et al 2013]() Equation 3
MCharnov2013 <- function(){
  self <- extend(Node,'MCharnov2013')

  self$predictand <- 'm'
  self$predictors <- c('linf','lmat','k')
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,k*((lmat/linf)^-1.5))
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
  
  self$tests <- function(){
  }
  
  self
}
