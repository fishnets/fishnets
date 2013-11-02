require(plyr)

#' Taxon.lookuppers are network nodes that are able to impute
#' missing taxonomic levels using simple table lookups. They are
#' used in fishnets to save the user having to do that themselves.
#' e.g. the user can just enter a species name
Taxon.lookupper <- function(from,to){
  self <- object('Taxon.lookupper')
  
  self$from <- from
  self$to <- to
  
  self$fit <- function(data){
    self$table <- ddply(data,self$from,function(sub){
        head(as.character(sub[,self$to]),n=1)
    })
    names(self$table) = c(self$from,self$to)
  }
  
  self$predict <- function(variables){
    which <- match(variables[[self$from]],self$table[,self$from])
    if(!sum(!is.na(which))==1) stop(paste("Key value is absent:",variables[[self$from]]))
    self$table[which,self$to]
  }
  
  self
}

Class.lookupper = function() Taxon.lookupper('order','class')
Order.lookupper = function() Taxon.lookupper('family','order')
Family.lookupper = function() Taxon.lookupper('genus','family')
Genus.lookupper = function() Taxon.lookupper('species','genus')
