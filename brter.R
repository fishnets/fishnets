require(gbm)
require(dismo)

#' A network node that predicts using Boosted Regression Trees (BRT)
#' 
#' @author Nokome Bentley
#' 
#' @param formula A formula for the terms to be included in the BRT
#' @param transform A function to apply to predicted value e.g. if formula is log(k)~... use exp
#' @param autotrees Should the number of trees be optimised usind gbm.step?
Brter <- function(formula,transform=identity,ntrees=2000,tree.complexity=10,learning.rate=0.001,
                  bag.fraction = 0.5,max.trees=5000){
  self <- extend(Node,'Brter')
  
  self$formula <- formula
  self$transform <- transform
  self$predictand <- all.vars(formula)[1]
  self$predictors <- all.vars(formula)[-1]
    
  self$fit <- function(data){
    # Create a model frame using the formula. This creates a
    # matrix that can be used for fitting the trees
    frame <- model.frame(self$formula,data)
    
    if(ntrees>0){
      self$brt <- gbm.fixed(
        data = frame,
        gbm.y = 1,
        gbm.x = 2:ncol(frame), 
        family = "gaussian",
        tree.complexity = tree.complexity,
        learning.rate = learning.rate,
        bag.fraction = bag.fraction,
        n.trees = ntrees
      )
      self$ntrees <- ntrees
    } else {
      self$brt <- gbm.step(
        data = frame,
        gbm.y = 1,
        gbm.x = 2:ncol(frame), 
        family = "gaussian",
        tree.complexity = tree.complexity,
        learning.rate = learning.rate,
        bag.fraction = bag.fraction,
        max.trees = max.trees
      )
      self$ntrees <- length(self$brt$trees)
    }
  }
  
  self$summary <- function(){
    print(summary(self$brt))
    #gbm.plot(self$brt)
    #gbm.plot.fits(self$brt)
  }
  
  # 'safe' predict function that only predicts
  # if predictand value is absent
  self$predict.safe <- function(data,transform=T,na.strict=T,na.keep=F) {
    data <- as.data.frame(data)
    
    # include predictand in data if it is missing
    if(!(self$predictand %in% names(data)))
      data[,self$predictand] <- NA
    
    # record locations that should not be
    # overwritten
    safe.loc <- !is.na(data[,self$predictand])
    
    # create a model frame using the formula
    # and keeping all NA values
    frame <- model.frame(self$formula,data,na.action=na.pass)
    
    # Generate predictions of length equal
    # to nrows(data)
    preds <- frame[,1]
    preds[!safe.loc] <- predict.gbm(
      self$brt,
      newdata = frame,
      n.trees = self$brt$gbm.call$best.trees,
      type = "response"
    )[!safe.loc]
    
    # Record s.d. of residuals for $sample()
    self$error <- sd(self$brt$residuals,na.rm=T)
    
    # if na.strict we suppose that we are not justified in predicting
    # predictand if any of the covariates are missing
    if(na.strict) {
      na.loc <- apply(frame,1,function(x) any(is.na(x[-1])))
      preds[!safe.loc &na.loc] <- NA
    }
    
    # if !na.keep remove all NA's from prediction vector
    if(!na.keep) preds <- preds[!is.na(preds)]
    
    # return with or without back-transformation
    # of predicted values
    if(transform) return(self$transform(preds))
    else return(preds)
  }
  
  self$predict <- function(data,transform=T,na.strict=T,na.keep=F){
    data <- as.data.frame(data)
     
    # include predictand in data if it is missing
    if(!(self$predictand %in% names(data)))
      data[,self$predictand] <- NA
    
    # create a model frame using the formula
    # and keeping all NA values
    frame <- model.frame(self$formula,data,na.action=na.pass)
    
    # Generate predictions of length equal
    # to nrows(data)
    preds <- predict.gbm(
      self$brt,
      newdata = frame,
      n.trees = self$brt$gbm.call$best.trees,
      type = "response"
    )
    # Record s.d. of residuals for $sample()
    self$error <- sd(self$brt$residuals,na.rm=T)
    
    # if na.strict we suppose that we are not justified in predicting
    # predictand if any of the covariates are missing
    if(na.strict) {
      na.loc        <- apply(frame,1,function(x) any(is.na(x[-1])))
      preds[na.loc] <- NA
    }
    
    # if !na.keep remove all NA's from prediction vector
    if(!na.keep) preds <- preds[!is.na(preds)]
    
    # return with or without back-transformation
    # of predicted values
    if(transform) return(self$transform(preds))
    else return(preds)
  }
  
  self$sample <- function(data){
    self$transform(rnorm(nrow(data),mean=self$predict(data,transform=F),sd=self$error))
  }
  
  self
}
