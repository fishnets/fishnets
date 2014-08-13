#' A `Node` for matural mortality based on
#' [Charnov et al 2013]() Equation 3
MCharnov2013Fitted <- function(){
  self <- extend(MCharnov2013,'MCharnov2013Fitted')

  formula <- log(m/k) ~ -1 + log(lmat/linf)
    
  self$fit <- function(data, ...){
    self$glm <- glm(formula,data=data)
  }
  
  self$predict <- function(data){
    logmk <- predict.glm(self$glm,newdata=data,type='response')
    exp(logmk)*data$k
  }
  
  self$n <- function(data) {
    # number of data points used for fitting
    frame <- model.frame(formula,data)
    nrow(frame)
  }
  
  self$sample <- function(data){
    # Get predictions with errors and no transformation
    predictions <- as.data.frame(predict.glm(self$glm,newdata=data,type='response',se.fit=T))
    # Calculate a standard deviation that combines se.fit and residual s.d.
    sigma <- sqrt(predictions$se.fit^2 + predictions$residual.scale^2)
    # Sample from normal distribution with that sigma
    preds <- rnorm(nrow(predictions),mean=predictions$fit,sigma)
    # Apply post transformation
    exp(preds)*data$k
  }
  
  self
}
