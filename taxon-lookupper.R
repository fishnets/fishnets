require(plyr)

#' Taxon.lookuppers are network nodes that are able to impute
#' missing taxonomic levels using simple table lookups. They are
#' used in fishnets to save the user having to do that themselves.
#' e.g. the user can just enter a species name
Taxon.lookupper <- function(from,to){
  self <- extend(Node,'Taxon.lookupper')
  
  self$predictors <- from
  self$predictand <- to
  
  self$fit <- function(data){
    self$table <- ddply(data,self$predictors,function(sub){
        head(as.character(sub[,self$predictand]),n=1)
    })
    names(self$table) = c(self$predictors,self$predictand)
  }
  
  self$predict <- function(data){
    which <- match(data[[self$predictors]],self$table[,self$predictors])
    if(!sum(!is.na(which))==1) stop(paste("Key value is absent:",data[[self$predictors]]))
    self$table[which,self$predictand]
  }
  
  self$sample <- function(data,samples=1){
    # Samples direct from prediction i.e. 100% certainty
    rep(self$predict(data),samples)
  }
  
  self
}

Class.lookupper = function() Taxon.lookupper('order','class')
Order.lookupper = function() Taxon.lookupper('family','order')
Family.lookupper = function() Taxon.lookupper('genus','family')
Genus.lookupper = function() Taxon.lookupper('species','genus')
