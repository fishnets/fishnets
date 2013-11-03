Species.any <- function(from,to){
  self <- extend(Node,'Species.any')
  
  self$predictors <- 'fish'
  self$predictand <- 'species'
  
  self$fit <- function(data){
    self$species <- unique(data$species)
  }
  
  self$predict <- function(data){
    sample(self$species,1)
  }
  
  self$sample <- function(data,samples=1){
    sample(self$species,samples)
  }
  
  self
}
