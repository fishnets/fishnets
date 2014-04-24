#' A network node that predicts based a Generalised Linear Model (GLM)
#' 
#' @author Nokome Bentley
#' 
#' @param formula The GLM formula
#' @param transform A function to apply to predicted value after predict.glm. e.g. if formula is log(k)~... use exp
#' @param factors.min Minimum records for a level of factor predictor variables
#' @param numerics.impute Logical flag indicating whether missing numeric predictans should be imputed
Glmer <- function(formula,transform=identity,factors.min=5){
  self <- extend(Node,'Glmer')
  
  self$formula <- formula
  self$transform <- transform
  self$factors.min <- factors.min
  
  self$predictand <- all.vars(self$formula)[1]
  self$predictors <- all.vars(self$formula)[-1]
  
  self$fit <- function(data){    
    # Restrict data to rows with variable of interest
    data <- data[!is.na(data[,self$predictand]),]
    # Prepare the data before fiting
    for(name in self$predictors){
      values <- data[,name]
      # Rationalise factors by collapsing poorly represented
      # levels into an '<other>' level and putting NAs into an '<unknown>' level
      if(is.factor(values) | is.character(values)){
        values <- as.character(values)
        # Calculate n in each level
        levels <- as.data.frame(table(values),responseName="n")
        # Determine main levels
        levels <- subset(levels,n>=self$factors.min)
        # Set others, unknowns and recreate factor
        values[!is.na(values) & !(values %in% levels$values)] <- '<other>'
        values[is.na(values)] <- '<unknown>'
        data[,name] <- factor(values,levels=unique(c(unique(values),'<other>','<unknown>')))
        # Can't fit model when there is less than 2 levels to a factor. In that case
        # drop it from the term
        uniques <- length(unique(data[,name]))
        if(uniques<2){
          self$formula <- update.formula(self$formula,paste(".~.-",name,sep=""))
          warning("Glmer removed term '",name,"' from formula because it had too few levels: ",uniques,". New formula: ",deparse(self$formula))
        }
      }
    }
    # Recalculate predictors in case any were removed
    self$predictors <- all.vars(self$formula)[-1]
    # Fit the GLM
    self$glm <- glm(self$formula,data=data)
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
        # If the orginal data did not have any records with a particular value for factor then prediction 
        # will not work (impossible for GLM to estimate a coefficient). So set those as NA to avoid an error
        counts <- table(self$glm$data[,name])
        values[!(values %in% names(counts[counts>0]))] <- NA
        values <- factor(values,levels=levels(self$glm$data[,name]))
        data[,name] <- values
      }
    }
    predictions <- predict.glm(self$glm,newdata=data,type='response',se.fit=errors)
    # When called with errors=T do not apply self$post transformation
    if(errors) return(as.data.frame(predictions))
    else return(self$transform(predictions))
  }
  
  #' Tune model formulas simple search over a range of values
  self$tune <-function(data,trials){
    results <- NULL
    for(formula in trials){
      cat("formula:",deparse(formula),"\n")
      self$formula = formula
      results <- rbind(results,data.frame(
        formula = deparse(formula),
        self$cross(data)$summary
      ))
    }
    list(
      best = results$formula[which.min(results$mpe)],
      trials = results
    )
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
