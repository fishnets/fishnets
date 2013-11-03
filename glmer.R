#' A network node that predicts based a Generalised Linear Model (GLM)
#' 
#' @param formula The GLM formula
#' @param post A function to apply to predicted value after predict.glm. e.g. if formula is log(k)~... use exp
#' @param numerics.min Minimum records for taxonomic imputation of numeric predictor variables
#' @param factors.min Minimum records for a level of factor predictor variables
Glmer <- function(formula,transform=identity,numerics.min=3,factors.min=5){
  self <- extend(Node,'Glmer')
  
  self$formula <- formula
  self$transform <- transform
  self$predictand <- all.vars(formula)[1]
  self$predictors <- all.vars(formula)[-1]
  self$numerics.min <- numerics.min
  self$factors.min <- factors.min
  
  self$fit <- function(data){
    # Restrict data to rows with variable of interest
    data <- data[!is.na(data[,self$predictand]),]
    # Rationalise factors by collapsing poorly represented
    # levels into an '<other>' level and putting NAs into an '<unknown>' level
    for(name in self$predictors){
      values <- data[,name]
      if(is.factor(values)){
        values <- as.character(values)
        # Calculate n in each level
        levels <- as.data.frame(table(values),responseName="n")
        # Determine main levels
        levels <- subset(levels,n>=self$factors.min)
        # Set others, unknowns and recreate factor
        values[!is.na(values) & !(values %in% levels$values)] <- '<other>'
        values[is.na(values)] <- '<unknown>'
        data[,name] = factor(values,levels=unique(c(unique(values),'<other>','<unknown>')))
      }
    }
    # Fit the GLM
    self$glm <- glm(formula,data=data)
  }
  
  self$predict <- function(data,errors=F){
    data <- as.data.frame(data)
    # Collapse factors down to the same levels as in the fitted GLM
    for(name in self$predictors){
      values <- data[,name]
      if(is.factor(values) | is.character(values)){
        values <- as.character(values)
        values[is.na(values)] <- '<unknown>'
        values[is.na(match(values,levels(self$glm$data[,name])))] <- '<other>'
        values <- factor(values,levels=levels(self$glm$data[,name]))
        data[,name] <- values
      }
    }
    predictions <- predict.glm(self$glm,newdata=data,type='response',se.fit=errors)
    # When called with errors=T do not apply self$post transformation
    if(errors) return(as.data.frame(predictions))
    else return(self$transform(predictions))
  }
  
  self$sample <- function(data,samples=1){
    # Get predictions with errors
    predictions <- self$predict(data,errors=T)
    # Calculate a standard deviation that combines se.fit and residual s.d.
    sigma <- sqrt(predictions$se.fit^2 + predictions$residual.scale^2)
    # Sample from normal distribution with that sigma
    preds <- rnorm(samples,predictions$fit,sigma)
    # Apply post transformation
    self$transform(preds)
  }
  
  self
}
