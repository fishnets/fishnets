Node <- function(){
  
  self <- object('Node')
  
  #' Cross validation of a node
  #' 
  #' @param data The data.frame to cross validate against
  #' @param folds Number of folds (number of times that data is split into training and testing datasets)
  #' @param fit Function for fitting model. Should return an object with a "predict" method.
  #' @param ... Other arguments to fit function
  self$cross <- function(folds,data){
    # Begin by removing all data rows with NAs for variable of interest
    data = data[!is.na(data[,self$predictand]),]
    # Randomly assign rows of data to a fold
    # Do this in a way that divides the data as evenly as possible..
    # Systematically assign folder number to give very close to even numbers in each fold...
    folder = rep(1:folds,length.out=nrow(data))
    # Randomly rearrange
    folder = sample(folder)
    # Result data.frame
    results = NULL
    # For each fold...
    for(fold in 1:folds){
      cat("Fold",fold,":")
      # Define training a testing datasets
      train = data[folder!=fold,]
      test = data[folder==fold,]
      # Fit model to training data
      self$fit(train)
      # Get predictions for testing data
      preds = self$predict(test)
      # Get the dependent variable
      tests = test[,self$predictand]
      # Calculate various prediction errors
      me = mean(abs(tests-preds))
      mse = mean((tests-preds)^2)
      mpe = mean(abs((tests-preds)/tests))
      r2 = cor(tests,preds)^2
      # Add to results
      results = rbind(results,data.frame(
        fold = fold,
        me = me,
        mse = mse,
        mpe = mpe,
        r2 = r2
      ))
      cat(mpe,r2,"\n")
    }
    
    # Summarise results
    summary = colMeans(results)
    summary = as.list(summary)
    summary$fold <- NULL
    
    return (list(
      summary = summary,
      folds = results
    ))
  }
   
   self
}




