#' A `Node` for matural mortality based on
#' [Pauly 1980]() Equation 3
MPauly1980Fitted <- function(){
  self <- extend(MPauly1980,'MPauly1980Fitted')

  formula <- log(m) ~ log(linf) + log(k) + log(temp)
    
  self$fit <- function(data, ...){
    self$glm <- glm(formula,data=data)
  }
  
  self$predict <- function(data){
    logm <- predict.glm(self$glm,newdata=data,type='response')
    exp(logm)
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
    preds <- suppressWarnings(rnorm(nrow(predictions),mean=predictions$fit,sigma))
    # Apply post transformation
    exp(preds)
  }
  
  self
}
