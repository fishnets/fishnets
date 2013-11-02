#'A network of predictor nodes
#'
#'@name Fishnet
#'@param nodes A list of network nodes
Fishnet <- function(
  nodes = list()
){
  self <- object('Fishnet')
  
  self$nodes <- nodes
  
  #'Fit the network using data
  #'
  #'@name Fishnet$fit
  #'@param data Data to be fitted to
  self$fit <- function(data){
    for(name in names(self$nodes)) self$nodes[[name]]$fit(data)
  }
  
  #'Predict values for variables from the network
  #'
  #'@name Fishnet$predict
  #'@param data A list of variables
  self$predict <- function(...){
    data = as.list(match.call())[-1]
    for(name in names(self$nodes)){
      if(is.null(data[[name]])){
        data[[name]] <- self$nodes[[name]]$predict(as.data.frame(data))
      }
    }
    data
  }
  
  self
}
