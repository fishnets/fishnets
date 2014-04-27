require(plyr)

#' A network node that looks up missing taxonomic levels
#' 
#' TaxonLookuppers are network nodes that are able to impute
#' missing taxonomic levels using simple table lookups. They are
#' used in Fishnets to save the user having to do that themselves.
#' e.g. the user can just enter a species name
#' 
#' @author Nokome Bentley
#' 
#' @param from Taxonomic level that is the key for lookup
#' @param to Taxonomic level that is the lookup value
TaxonLookupper <- function(from,to){
  self <- extend(Node,'TaxonLookupper')
  
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
  
  self$sample <- function(data){
    # Samples direct from prediction i.e. 100% certainty
    self$predict(data)
  }
  
  self
}

ClassLookupper = function() extend(TaxonLookupper,'ClassLookupper','order','class')
OrderLookupper = function() extend(TaxonLookupper,'OrderLookupper','family','order')
FamilyLookupper = function() extend(TaxonLookupper,'FamilyLookupper','genus','family')
GenusLookupper = function() extend(TaxonLookupper,'GenusLookupper','species','genus')
