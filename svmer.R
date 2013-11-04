require(e1071)

#' A network node that predicts based on a Support Vector Machine (SVM)
#' 
#' @author Nokome Bentley
#' 
#' @param formula The SVM formula
#' @param transform A function to apply to predicted value e.g. if formula is log(k)~... use exp
Svmer <- function(formula,transform=identity){
  self <- extend(Node,'Svmer')
  
  self$formula <- formula
  self$transform <- transform
  self$predictand <- all.vars(formula)[1]
  self$predictors <- all.vars(formula)[-1]
  
  self$fit <- function(data){
    # Record levels of factors so that the same levels can be applied in $predict()
    self$levels <- list()
    for(name in self$predictors){
      if(is.factor(data[[name]])) self$levels[[name]] <- levels(data[[name]])
    }
    # Fit da model
    self$svm <- svm(formula,data=data,cost=1000,gamma=0.0001)
    # Record s.d. of residuals for $sample()
    self$error <- sd(self$svm$residuals)
  }
  
  self$predict <- function(data,transform=T){
    data <- as.data.frame(data)
    # Apply levels of factors used in fitting
    for(name in self$predictors){
      if(is.factor(data[[name]]) | is.character(data[[name]])) data[[name]] <- factor(data[[name]],levels=self$levels[[name]])
    }
    # Get da predictions
    preds <- predict(self$svm,newdata=data[,self$predictors])
    if(transform) return(self$transform(preds))
    else return(preds)
  }
  
  self$sample <- function(data,samples=1){
    self$transform(rnorm(samples,self$predict(data,transform=F),self$error))
  }
  
  self
}
