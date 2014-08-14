#' A `Node` for matural mortality based on
#' [Charnov et al 2013]() Equation 3
MHoenig1983Fitted <- function(){
  self <- extend(MHoenig1983,'MHoenig1983Fitted')

  formula <- log(m) ~ log(amax)
  
  self$fit <- function(data, ...){
    self$glm <- glm(formula,data=data)
  }
  
  self$predict <- function(data,transform=T,na.strict=T,na.keep=T){
    
    data <- as.data.frame(data)
    
    # by default predict.glm() will predict NA for
    # any data row with missing covariate values
    # consistent with na.strict=T
    preds <- predict.glm(self$glm,newdata=data,type='response')
    preds <- exp(preds)
    
    # if !na.keep remove all NA's from prediction vector
    if(!na.keep) preds <- preds[!is.na(preds)]
    
    return(preds)
  }
  
  self$predict.safe <- function(data,transform=T,na.strict=T,na.keep=T) {
    
    data <- as.data.frame(data)
    
    if(self$predictand %in% names(data)) {
      safe.loc <- !is.na(data[,self$predictand])
    } else {
      safe.loc <- !numeric(nrow(data))
    }
    
    # by default predict.glm() will predict NA for
    # any data row with missing covariate values
    # consistent with na.strict=T
    preds <- predict.glm(self$glm,newdata=data,type='response')
    preds <- exp(preds)
    
    # restore existent values
    preds[safe.loc] <- data[safe.loc,self$predictand]
    
    # if !na.keep remove all NA's from prediction vector
    if(!na.keep) preds <- preds[!is.na(preds)]
    
    return(preds)
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
