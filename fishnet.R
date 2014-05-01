#' A network of predictor nodes
#'
#' @author Nokome Bentley
#'
#' @param nodes A list of network nodes
Fishnet <- function(...){
  self <- object('Fishnet')
  
  # Convert call arguments to a list of nodes
  self$nodes <- list(...)
  
  #' Fit the network using data
  #'
  #' @name Fishnet$fit
  #' @param data Data to be fitted to
  self$fit <- function(data){
    # Delete stored results since we are refitting network
    self$stored <- NULL
    # Fit all nodes
    for(name in names(self$nodes)){
      cat('Fitting',name)
      self$nodes[[name]]$fit(data)
      cat('.\n')
    }
    self
  }
  
  #' Predict values for variables from the network
  #'
  #' @name Fishnet$predict
  #' @param data A list of variable valuables
  self$predict <- function(data){
    results <- as.data.frame(data)
    # Get nodes to predict values for their predictands
    for(name in names(self$nodes)){
      if(is.null(results[[name]])){
        results[,name] <- self$nodes[[name]]$predict(results)
      }
    }
    results
  }
  
  #' Sample the network based on a list or data.frame of variable values
  #' 
  #' @name Fishnet$df
  #' @param data A data.frame or list of variable values
  #' @param samples Number of samples to generate
  self$sample_df <- function(data,samples=1){
    data <- as.data.frame(data)
    # Problems occur if there is only one column in the df
    # so add a dummy one that gets taken out later
    data$dummy <- NA
    # Repeat the data.frame `samples` times so that the `sample` method
    # of each node can operate in a vectorised way
    results <- data[rep(1:nrow(data),samples),]
    # Iterate over nodes filling in values as necessary
    for(name in names(self$nodes)){
      if(!(name %in% names(results))){
        results[,name] <- self$nodes[[name]]$sample(results)
      }
    }
    # Reove that dumny column
    results$dummy <- NULL
    results
  }
  
  #' Sample the network based on a list of distributions
  #' 
  #' @name Fishnet$sample_distributions
  #' @param dists A list of Distributions
  #' @param samples Number of samples to generate
  self$sample_distributions <- function(dists,samples=1){
    # Create a data frame with the right number of rows
    results <- data.frame(dummy=rep(NA,samples))
    # Add random samples from each distribution
    for(name in names(dists)) results[,name] <- dists[[name]]$random(nrow(results))
    results$dummy <- NULL
    # Iterate over nodes filling in values as necessary
    for(name in names(self$nodes)){
      if(!(name %in% names(results))){
        results[,name] <- self$nodes[[name]]$sample(results)
      }
    }
    results
  }  
  
  #' Sample the network
  #' 
  #' @name Fishnet$sample
  #' @param from A list of values or a list of Distributions to sample from
  #' @param samples Number of samples to generate
  self$sample <- function(from,samples=10000){
    if(inherits(from,'Distributions')){
      return (self$sample_distributions(from,samples))
    }
    if(inherits(from,'data.frame') | inherits(from,'list')){
      return (self$sample_df(from,samples))
    }
    stop(paste('Unable to handle data of type:',paste(class(data),collapse=",")))
  }
  
  #' Create a graph visualization of the network
  #' Requires graphviz to be installed
  #' 
  #' @name Fishnet$graph
  #' @param folder Directory in which graph will be produced
  self$graph <- function(folder="."){
    # Create the folder and DOT file
    dir.create(folder,showWarnings=F,recursive=T)
    dotname <- file.path(folder,"graph.dot")

    out <- file(dotname,open='w')
    write <- function(...) cat(...,file=out)
    write('
    strict digraph G {
      node [shape=box,fontsize=9,fontname=Helvetica];
      edge [fontsize=9,fontname=Helvetica];
    ')
    
    for(name in names(self$nodes)){
      node <- self$nodes[[name]]
      write(node$predictand,';\n')
      for(predictor in node$predictors){
        write(predictor,'->',node$predictand,';\n')
      }
    }
    
    write('}')
    close(out)
    
    system(paste('dot -Tsvg ',dotname,' -o ',file.path(folder,"graph.svg")))
    system(paste('dot -Tpng ',dotname,' -o ',file.path(folder,"graph.png")))
  }
  
  #' Store samples from this Fishnet
  #' 
  #' @param from Data to be used as starting point of sampling
  #' @param samples Number of samples
  #' @param folder Directory where the samples will be stored
  #' 
  #' @note To be implemented: as well as storing to .RData this method should store to a text or binary
  #' file that is suitable for reading into C++ or other language
  self$store <- function(from,samples,folder){
    # Generate samples
    samples <- self$sample(from,samples)
    # Create the folder
    dir.create(folder,showWarnings=F,recursive=T)
    # Save the samples
    save(samples,file=file.path(folder,"store.RData"))
    # Save documentation on this Fishnet
    self$graph(folder)
  }
  
  #' Restore samples into this Fishnet
  #' 
  #' Rather than fitting a Fishnet it can be restored from a previous call to $store()
  #' 
  #' @param folder Directory where the samples will be stored
  #' 
  #' @note To be implemented: methods $sample() and $predict() should check for is.null(self$stored) and sample
  #' from stored if present.
  self$restore <- function(folder){
    # Load the store and assign to self
    load(file.path(folder,"store.RData"))
    self$stored = samples
  }
  
  self
}
