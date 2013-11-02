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
  #'@param variables A list of variables
  self$predict <- function(variables){
    if(!is.list(variables)) stop('Argument "variables" must be a list')
    for(name in names(self$nodes)){
      if(is.null(variables[[name]])) variables[[name]] <- self$nodes[[name]]$predict(variables)
    }
    variables
  }
  
  self
}
