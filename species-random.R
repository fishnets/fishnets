#' A simple network node that 'predicts' the species.
#' This node simply draws a random species names from the fitted data.
#' It is needed for the network to produce a generic prior (i.e not specific to any species)
#'
#' @author Nokome Bentley
SpeciesRandom <- function(){
  self <- extend(Node,'SpeciesRandom')
  
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
