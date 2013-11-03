#' A network of predictor nodes
#'
#' @name Fishnet
#' @param nodes A list of network nodes
Fishnet <- function(
  nodes = list()
){
  self <- object('Fishnet')
  
  self$nodes <- nodes
  
  #' Fit the network using data
  #'
  #' @name Fishnet$fit
  #' @param data Data to be fitted to
  self$fit <- function(data){
    for(name in names(self$nodes)) self$nodes[[name]]$fit(data)
  }
  
  #' Predict values for variables from the network
  #'
  #' @name Fishnet$predict
  #' @param data A list of variable valuables
  self$predict <- function(data){
    # Get nodes to predict values for their predictands
    for(name in names(self$nodes)){
      if(is.null(data[[name]])){
        data[[name]] <- self$nodes[[name]]$predict(as.data.frame(data))
      }
    }
    data
  }
  
  #' Sample the network based on a list of variable values
  #' 
  #' @name Fishnet$sample.list
  #' @param data A list of variable valuables
  #' @param samples Number of samples to generate
  self$sample.list <- function(data,samples=1){
    # Get nodes to provide one sample and then pass this on to remaining
    # nodes thus propogating prediction errors through the network
    results <- NULL
    for(sample in 1:samples){
      sample_data <- as.data.frame(data)
      for(name in names(self$nodes)){
        if(is.null(sample_data[[name]])){
          sample_data[[name]] <- self$nodes[[name]]$sample(sample_data,1)
        }
      }
      results <- rbind(results,sample_data)
    }
    results
  }
  
  self$sample.distributions <- function(dists,samples=1){
    results <- NULL
    for(sample in 1:samples){
      sample_data <- as.data.frame(lapply(dists,function(dist) dist$random()))
      for(name in names(self$nodes)){
        if(is.null(sample_data[[name]])){
          sample_data[[name]] <- self$nodes[[name]]$sample(sample_data,1)
        }
      }
      results <- rbind(results,sample_data)
    }
    results
  }  
  
  self$sample <- function(arg,samples=1){
    if(inherits(arg,'Distributions')){
      return (self$sample.distributions(arg,samples))
    }
    if(inherits(arg,'list')){
      return (self$sample.list(arg,samples))
    }
    stop(paste('Unable to handle data of type:',paste(class(data),collapse=",")))
  }
  
  self
}
