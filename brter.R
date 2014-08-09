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
  
  self$formula    <- formula
  self$transform  <- transform
  self$predictand <- all.vars(formula)[1]
  self$predictors <- all.vars(formula)[-1]
  
  if(length(self$predictors)<2) 
    stop('Formula must contain at least 2 predictors\n')
      
  self$fit <- function(data,pars){
    
    # Create a model frame using the formula. This creates a
    # matrix that can be used for fitting the trees
    frame <- suppressWarnings(model.frame(self$formula,as.data.frame(data)))
    
    # make sure character vectors in frame are factors
    for(par in names(frame)) 
      if(is.character(frame[,par]))
        frame[,par] <- as.factor(frame[,par])
    
    # the 'pars' argument allows further control of which
    # covariates should be included from the formula. It defaults
    # to 'all covariates' and can be specified as either a character 
    # vector or numeric reference to particular columns. This is useful
    # for selection of predictors
    if(missing(pars)) { pars <- 2:ncol(frame)
    } else {
      if(is.character(pars)) pars <- match(pars,names(frame))
    }
    
    if(length(pars)<2) 
      stop('Formula must contain at least 2 predictors\n')
    
    if(ntrees>0){
      self$brt <- gbm.fixed(
        data = frame,
        gbm.y = 1,
        gbm.x = pars, 
        family = "gaussian",
        tree.complexity = tree.complexity,
        learning.rate = learning.rate,
        bag.fraction = bag.fraction,
        n.trees = ntrees
      )
      self$ntrees <- ntrees
      self$n <- length(self$brt$fit)
    } else {
      self$brt <- gbm.step(
        data = frame,
        gbm.y = 1,
        gbm.x = pars, 
        family = "gaussian",
        tree.complexity = tree.complexity,
        learning.rate = learning.rate,
        bag.fraction = bag.fraction,
        max.trees = max.trees
      )
      self$ntrees <- length(self$brt$trees)
      self$n <- length(self$brt$fit)
    }
    cat('predictors:',self$brt$gbm.call$predictor.names,'\n')
  }
  
  self$summary <- function(){
    print(summary(self$brt))
    #gbm.plot(self$brt)
    #gbm.plot.fits(self$brt)
  }
  
  # 'safe' predict function that only predicts
  # if data value is absent
  self$predict.safe <- function(data,transform=T,na.strict=T,na.keep=F) {
    data <- as.data.frame(data)
    # Create a model frame using the formula. This creates a
    # data.frame that can be used for fitting the trees
    
    if(!(self$predictand %in% names(data)))
      data[,self$predictand] <- NA
    safe.loc <- !is.na(data[,self$predictand])
            
    # create data frame
    # for model prediction
    frame <- model.frame(self$formula,data,na.action=na.pass)

    # Generate predictions from frame
    # using predictors recorded in gbm.object
    preds <- frame[,1]
    preds[!safe.loc] <- predict.gbm(
                                  self$brt,
                                  newdata = frame,
                                  n.trees = self$brt$gbm.call$best.trees,
                                  type = "response"
                                )[!safe.loc]
    # Record s.d. of residuals for $sample()
    self$error <- sd(self$brt$residuals)
    
    # if na.strict we suppose that we are not justified in predicting
    # predictand if any of the covariates is NA
    if(na.strict) {
      na.loc <- apply(frame,1,function(x) any(is.na(x[-1])))
      preds[!safe.loc &na.loc] <- NA
    }
    
    # if !na.keep remove all NA's from prediction vector
    if(!na.keep) preds <- preds[!is.na(preds)]
    
    if(transform) return(self$transform(preds))
    else return(preds)
  }
  
  self$predict <- function(data,transform=T,na.strict=T,na.keep=F){
    data <- as.data.frame(data)
    # Create a model frame using the formula. This creates a
    # data.frame that can be used for fitting the trees
    
    if(!(self$predictand %in% names(data)))
      data[,self$predictand] <- NA
    
    #if(na.keep & !nrow(data)>1) 
    #  na.keep <- FALSE
    #if(na.keep) {
    #  # if keeping NA's 
    #  # 1. create output vector of correct length
    #  # 2. populate with current values of predictand
    #  # 3. overwrite current values if data available
    #  
    #  # create vector of locations for which prediction will
    #  # not be preformed (na.loc) and prediction vector of correct
    #  # length with populated values (these will be overwritten 
    #  # at na.preds[!na.loc])
    #  frame.pass <- model.frame(self$formula,data,na.action=na.pass)
    #  na.loc     <- apply(frame.pass,1,function(x) any(is.na(x[-1])))
    #  na.preds   <- frame.pass[,1]
    #}
    
    # create data frame with all NA's removed
    # for model prediction
    #data[,self$predictand] <- 1
    frame <- model.frame(self$formula,data,na.action=na.pass)
    # Generate predictions from frame
    # using predictors recorded in gbm.object (which
    # may be a subset of the data columns if 'pars' 
    # argument was included during fitting)
    preds <- predict.gbm(
      self$brt,
      newdata = frame,
      n.trees = self$brt$gbm.call$best.trees,
      type = "response"
    )
    # Record s.d. of residuals for $sample()
    self$error <- sd(self$brt$residuals)
    
    # if na.strict we suppose that we are not justified in predicting
    # predictand if any of the covariates is missing
    if(na.strict) {
      na.loc        <- apply(frame,1,function(x) any(is.na(x[-1])))
      preds[na.loc] <- NA
    }
    
    # if !na.keep remove all NA's from prediction vector
    if(!na.keep) preds <- preds[!is.na(preds)]
    
    #if(na.keep) {
    #  na.preds[!na.loc] <- preds  
    #  preds <- na.preds
    #}
    
    if(transform) return(self$transform(preds))
    else return(preds)
  }
  
  self$sample <- function(data){
    self$transform(rnorm(nrow(data),mean=self$predict(data,transform=F),sd=self$error))
  }
  
  self
}
