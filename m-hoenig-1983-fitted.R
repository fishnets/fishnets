#' A `Node` for matural mortality based on
#' [Charnov et al 2013]() Equation 3
MHoenig1983Fitted <- function(){
  self <- extend(MHoenig1983,'MHoenig1983Fitted')

  self$fit <- function(data){
    self$lm <- lm(log(m) ~ log(amax),data=data)
  }
  
  self$predict <- function(data){
    logm <- predict.lm(self$lm,newdata=data,type='response')
    exp(logm)
  }
  
  self$n <- function(data) {
    # number of data points used for fitting
  }
  
  self$sample <- function(data){
    # Get predictions with errors and no transformation
    predictions <- as.data.frame(predict.lm(self$lm,newdata=data,type='response',se.fit=T))
    # Calculate a standard deviation that combines se.fit and residual s.d.
    sigma <- sqrt(predictions$se.fit^2 + predictions$residual.scale^2)
    # Sample from normal distribution with that sigma
    preds <- suppressWarnings(rnorm(nrow(predictions),mean=predictions$fit,sigma))
    # Apply post transformation
    exp(preds)
  }
  
  self
}
