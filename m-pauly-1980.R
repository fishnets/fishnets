#' A `Node` for matural mortality based on
#' [Pauly 1980]()
MPauly1980 <- function(){
  self <- extend(Node,'MPauly1980')

  self$predictand <- 'm'
  self$predictors <- c('linf','k','temp')
  
  self$fit <- function(data, ...){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,exp(-0.0066-0.279*log(linf)+0.6543*log(k)+0.4634*log(temp)))
  }
  
  self$sample <- function(data){
    predictions <- log(self$predict(data))
    exp(rnorm(length(predictions),mean=predictions,sd=self$error))
  }
    
  self$tests <- function(){
  }
  
  self
}
