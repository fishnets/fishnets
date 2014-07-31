require(INLA)

#' A network node that predicts based a bayesian model estiamted using Laplace approximation (INLA)
#' 
#' @author Philipp Neubauer
#' 
#' @param formula The INLA formula
#' @param transform A function to apply to predicted value after predict.glm. e.g. if formula is log(k)~... use exp
#' @param factors.min Minimum records for a level of factor predictor variables
#' @param numerics.impute Logical flag indicating whether missing numeric predictans should be imputed
Bayser <- function(formula,transform=identity,factors.min=5){
  self <- extend(Node,'Bayser')
  
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
          warning("Bayser removed term '",name,"' from formula because it had too few levels: ",uniques,". New formula: ",deparse(self$formula))
        }
      }
    }
    # Recalculate predictors in case any were removed
    self$predictors <- all.vars(self$formula)[-1]
    # Fit the GLM
    self$inla <- inla(self$formula,data=data)
    self$fit_data <- data
  }
  
  self$predict <- function(data,errors=F,transform=T){
    
    data[,self$predictand] <- 1
    
    # subset the fit data to only model components; necessary to concatenate with prediction data
    fdata <- model.frame(paste(self$predictand,'~',paste(self$predictors,collapse='+')),self$fit_data)
    pdata <- model.frame(paste(self$predictand,'~',paste(self$predictors,collapse='+')),data)
    
    pdata[,self$predictand] <- NA
    
    ld <- nrow(pdata)
    data <- as.data.frame(rbind(pdata,fdata))
    for(name in self$predictors){
      values <- data[,name]
      if(is.factor(values) | is.character(values)){
        values <- as.character(values)
        values[is.na(values)] <- '<unknown>'
        values[is.na(match(values,levels(self$fit_data[,name])))] <- '<other>'
        values <- factor(values,levels=levels(self$fit_data[,name]))
        data[,name] <- values
      }
      uniques <- length(unique(data[,name]))
      if(uniques<2){
        self$formula <- update.formula(self$formula,paste(".~.-",name,sep=""))
        warning("Bayser removed term '",name,"' from formula because it had too few levels: ",uniques,". New formula: ",deparse(self$formula))
      }
    }
    predictions <- inla(self$formula,data=data,control.predictor=list(compute=T))$summary.linear.predictor[1:ld,1:2]
    if(errors) return(as.data.frame(predictions))
    if(transform) self$transform(predictions$mean)
    else predictions$mean
  }
  
  self$sample <- function(data){
    # Get predictions with errors and no transformation
    predictions <- self$predict(data,errors=T,transform=F)  
    # Sample from normal distribution with that sigma
    preds <- rnorm(nrow(predictions),mean=predictions$mean,predictions$sd)
    # Apply post transformation
    self$transform(preds)
  }
  
  self
}
