#' A `Node` for matural mortality based on
#' [Then et al 2014]() 
MThenEtAl2014Fitted <- function(){
  self <- extend(MThenEtAl2014,'MThenEtAl2014Fitted')

  self$fit <- function(data){
    self$glm <- glm(log(m*amax) ~ 1,data=data)
  }
  
  self$predict <- function(data){
    logmamax <- predict.glm(self$glm,newdata=data,type='response')
    exp(logmamax) / data$amax
  }
  
  self$n <- function(data) {
    # number of data points used for fitting
  }
  
  self$sample <- function(data){
    # Get predictions with errors and no transformation
    predictions <- as.data.frame(predict.glm(self$glm,newdata=data,type='response',se.fit=T))
    # Calculate a standard deviation that combines se.fit and residual s.d.
    sigma <- sqrt(predictions$se.fit^2 + predictions$residual.scale^2)
    # Sample from normal distribution with that sigma
    preds <- suppressWarnings(rnorm(nrow(predictions),mean=predictions$fit,sigma))
    # Apply post transformation
    exp(preds) / data$amax
  }
  
  self
}
