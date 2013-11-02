#' A network node that predicts based on taxonomic group.
#' It moves up the taxonomic levels until it finds a level which has at least nmin
#' values for the variable and then uses the median (for numeric variables) or mode (for factors etc)
#' at that taxonomic level.
#' 
#' @param variable The name of the variable of interest
#' @param nmin Minimum sample size required for median/mode
Taxonomic.imputer <- function(variable,nmin){
  self <- object('Taxonomic.imputer')
  
  self$variable <- variable
  self$nmin <- nmin
  
  self$fit <- function(data){
    # Determine whether to use median or mode for defualt based on type of variable...
    if(is.numeric(data[,self$variable])){
      defaults <- function(subset){    
        list(
          n = sum(!is.na(subset[,self$variable])),
          value = median(subset[,self$variable],na.rm=T)
        )
      }
    } else {
      defaults <- function(subset){
        list(
          n = sum(!is.na(subset[,self$variable])),
          value = names(which.max(table(subset[,self$variable])))
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
  
  self$predict <- function(variables){
    # Try each taxonomic level until sufficient sample size for the
    # variable is obtained
    for(level in c('species','genus','family','order','class')){
      value = as.character(variables[[level]])
      n <- tryCatch(self$defaults[[c(level,value,'n')]],error=function(error) 0)
      if(!is.null(n)){
        if(n>=self$nmin) return(self$defaults[[c(level,value,'value')]])
      }
    } 
    # If none of the taxonomic levels has sufficient sample size then just
    # return the default for all fish
    return(self$defaults$fish$all$value)
  }
  
  self
}
