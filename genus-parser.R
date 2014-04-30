# A `Node` to parse the genus name from a full
# binomial species name. Useful to use instead of
# a `SpeciesLookupper` node
GenusParser <- function(from,to){
  self <- extend(Node,'GenusParser')
  
  self$predictors <- 'species'
  self$predictand <- 'genus'
  
  self$predict <- function(data){
    unlist(lapply(strsplit(as.character(
      data$species
    ),' '),function(x)x[1]))
  }
  
  self$sample <- function(data){
    # Samples direct from prediction i.e. 100% certainty
    self$predict(data)
  }
  
  self
}
