
###############
# CV function #
###############

cvfishnet <- function(self,data,folds=10,byspecies=F) {
  
  if(byspecies) {
    data$species <- factor(data$species)
    spp <- levels(data$species)
    folds <- length(spp)
    levels(data$species) <- 1:folds
    folder <- as.numeric(data$species)
    levels(data$species) <- spp
  } else {
    folder = rep(1:folds,length.out=nrow(data))
    folder = sample(folder)
  }
  
  # Result data.frame
  results = NULL
  # For each fold...
  for(fold in 1:folds){
    #cat("Fold",fold,"\n")
    # Define training a testing datasets
    train = data[folder!=fold,]
    test = data[folder==fold,]
    
    # test vector
    tests <- test[,'m']
    test[,'m'] <- NA
    
    # Fit model to training data
    self$fit(train)
    
    # remove predictor columns with all NAs
    # in test data set (these will be imputed,
    # if some data are present no imputation will
    # be performed and predictand will be NA
    # for those rows with missing predictand)
    na.loc <- apply(test,2,function(x) all(is.na(x)))
    test   <- test[,!na.loc] 
    
    # impute missing predictor columns and
    # predict predictand
    for(name in names(self$nodes)){
      if(is.null(test[[name]])){
        test[,name] <- self$nodes[[name]]$predict(test)
      }
    }
    
    # get predictions
    preds <- test[,'m']
    
    # remove NA's
    na.loc <- is.na(preds) | is.na(tests)
    preds <- preds[!na.loc]
    tests <- tests[!na.loc]
    
    # Calculate various prediction errors
    if(!all(na.loc)) {
      me = mean(abs(tests-preds))
      mse = mean((tests-preds)^2)
      mpe = mean(abs((tests-preds)/tests))
      r2 = cor(tests,preds,use="pairwise.complete.obs")^2
      dev = mean((tests - preds) * (tests - preds))
      
      # Add to results
      results = rbind(results,data.frame(
        fold = fold,
        obs = mean(tests),
        hat = mean(preds),
        me = me,
        mse = mse,
        mpe = mpe,
        r2 = r2,
        dev = dev
      ))
    }
    
  }
  # Summarise results
  summary <- data.frame(mean=apply(results[,4:8],2,mean,na.rm=T),se=apply(results[,4:8],2,function(x) sd(x)/sqrt(length(x))))
  # Return summary and raw results
  return (list(
    summary = summary,
    folds = results
  ))
}
