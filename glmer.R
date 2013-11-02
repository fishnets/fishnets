#' A network node that predicts based a Generalised Linear Model (GLM)
#' 
#' @param formula The GLM formula
#' @param post A function to apply to predicted value after predict.glm. e.g. if formula is log(k)~... use exp
#' @param numerics.min Minimum records for taxonomic imputation of numeric predictor variables
#' @param factors.min Minimum records for a level of factor predictor variables
Glmer <- function(formula,post=identity,numerics.min=3,factors.min=5){
  self <- extend(Node,'Glmer')
  
  self$formula <- formula
  self$post <- post
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
  
  self$predict <- function(data){
    data <- as.data.frame(data)
    # Collapse factors down to the same levels as in the fitted GLM
    for(name in self$predictors){
      values <- data[,name]
      if(is.factor(values) | is.character(values)){
        values <- as.character(values)
        values[is.na(values)] <- '<unknown>'
        values <- factor(values,levels=levels(self$glm$data[,name]))
        values[is.na(values)] <- '<other>'
        data[,name] <- values
      }
    }
    self$post(predict.glm(self$glm,newdata=data))
  }
  
  self
}
