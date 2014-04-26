require(gbm)
require(dismo)

#' A network node that predicts using Boosted Regression Trees (BRT)
#' 
#' @author Nokome Bentley
#' 
#' @param formula A formula for the terms to be included in the BRT
#' @param transform A function to apply to predicted value e.g. if formula is log(k)~... use exp
#' @param autotrees Should the number of trees be optimised usind gbm.step?
Brter <- function(formula,transform=identity,ntrees=2000){
  self <- extend(Node,'Brter')
  
  self$formula <- formula
  self$transform <- transform
  self$predictand <- all.vars(formula)[1]
  self$predictors <- all.vars(formula)[-1]
  self$ntrees <- ntrees
    
  self$fit <- function(data){
    # Create a model frame using the formula. This creates a
    # matrix that can be used for fitting the trees
    frame <- model.frame(self$formula,data)
    
    if(self$ntrees>0){
      self$brt <- gbm.fixed(
        data = frame,
        gbm.y = 1,
        gbm.x = 2:ncol(frame), 
        family = "gaussian",
        tree.complexity = 10,
        learning.rate = 0.001,
        bag.fraction = 0.5,
        n.trees = self$ntrees
      )
    } else {
      self$brt <- gbm.step(
        data = frame,
        gbm.y = 1,
        gbm.x = 2:ncol(frame), 
        family = "gaussian",
        tree.complexity = 10,
        learning.rate = 0.001,
        bag.fraction = 0.5,
        max.trees = 5000
      )
    }
  }
  
  self$summary <- function(){
    print(summary(self$brt))
    #gbm.plot(self$brt)
    gbm.plot.fits(self$brt)
  }
  
  self$predict <- function(data,transform=T){
    data <- as.data.frame(data)
    # Create a model frame using the formula. This creates a
    # matrix that can be used for fitting the trees
    frame <- model.frame(self$formula,data)
    # Generate predictions
    preds <- predict.gbm(
      self$brt,
      frame,
      n.trees = self$brt$gbm.call$best.trees,
      type = "response"
    )
    if(transform) return(self$transform(preds))
    else return(preds)
  }
  
  self$sample <- function(data,samples=1){
    #self$transform(rnorm(samples,self$predict(data,transform=F),self$error))
  }
  
  self
}
