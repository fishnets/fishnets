#' A network node that predicts based on taxonomic group.
#' It moves up the taxonomic levels until it finds a level which has at least nmin
#' values for the variable and then uses the median (for numeric variables) or mode (for factors etc)
#' at that taxonomic level.
#' 
#' @param predictand The name of the variable of being predicted
#' @param nmin Minimum sample size required for median/mode
Taxonomic.imputer <- function(predictand,nmin){
  self <- extend(Node,'Taxonomic.imputer')
  
  self$predictand <- predictand
  self$predictors <- c('species','genus','family','order','class')
  self$nmin <- nmin
  
  self$fit <- function(data){
    # Determine whether to use median or mode for defualt based on type of predictand
    if(is.numeric(data[,self$predictand])){
      defaults <- function(subset){    
        list(
          n = sum(!is.na(subset[,self$predictand])),
          value = median(subset[,self$predictand],na.rm=T)
        )
      }
    } else {
      defaults <- function(subset){
        list(
          n = sum(!is.na(subset[,self$predictand])),
          value = names(which.max(table(subset[,self$predictand])))
        )
      }
    }
    # Calculate defaults for all taxomomic levels including all fish...
    self$defaults = list(
      fish = by(data,list(fish=rep("all",nrow(data))),defaults),
      class = by(data,list(class=data$class),defaults),
      order = by(data,list(order=data$order),defaults),
      family = by(data,list(family=data$family),defaults),
      genus = by(data,list(genus=data$genus),defaults),
      species = by(data,list(species=data$species),defaults)
    )
  }
  
  self$predict <- function(data){
    # Try each taxonomic level until sufficient sample size for the
    # variable is obtained
    preds <- NULL
    for(row in 1:nrow(data)){
      for(level in c('species','genus','family','order','class','fish')){
        group = as.character(data[row,level])
        n <- tryCatch(self$defaults[[c(level,group,'n')]],error=function(error) 0)
        if(!is.null(n)){
          if(n>=self$nmin){
            preds = c(preds,self$defaults[[c(level,group,'value')]][1])
            break
          }
        }
      } 
    }
    preds
  }
  
  self
}
