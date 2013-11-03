#' A network node that predicts based on taxonomic group.
#' It moves up the taxonomic levels until it finds a level which has at least nmin
#' values for the variable and then uses the median (for numeric variables) or mode (for factors etc)
#' at that taxonomic level.
#' 
#' @param predictand The name of the variable of being predicted
#' @param transform A pair of from/to transformation functions e.g. c(log,exp)
#' @param nmin Minimum sample size required for median/mode
Taxonomic.imputer <- function(predictand,transform,nmin){
  self <- extend(Node,'Taxonomic.imputer')
  
  self$predictand <- predictand
  self$transform <- transform
  self$predictors <- c('species','genus','family','order','class')
  self$nmin <- nmin
  
  self$fit <- function(data){
    # Determine the type of the predictand
    self$type <- if(is.numeric(data[,self$predictand])) numeric else factor
    # Determine type of default values to record based on type of predictand
    if(is.numeric(self$type())){
      defaults <- function(subset){  
        values <- self$transform[[1]](subset[,self$predictand])
        list(
          n = sum(!is.na(values)),
          expect = mean(values,na.rm=T),
          sd = sd(values,na.rm=T)
        )
      }
    } else {
      defaults <- function(subset){
        values <- subset[,self$predictand]
        list(
          n = sum(!is.na(values)),
          expect = names(which.max(table(values)))
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
  
  self$predict <- function(data,errors=F){
    preds <- NULL
    # Loop through data rows
    for(row in 1:nrow(data)){
      # Try each taxonomic level until sufficient sample size for the
      # variable is obtained
      best <- NULL
      for(level in c('species','genus','family','order','class','fish')){
        group <- as.character(data[row,level])
        defaults <- self$defaults[[c(level,group)]]
        if(!is.null(defaults)){
          if(defaults$n>=self$nmin){
            best <- defaults
            break
          }
        }
      }
      if(errors) return(best)
      else preds = c(preds,self$transform[[2]](best$expect))
    }
    preds
  }
  
  self$sample <- function(data,samples=1){
    # Currently, only provides sampling for numeric predictands
    if(!is.numeric(self$type())) stop('Sorry, sampling of non-numeric predictands is not implemented yet')
    # Get prediction
    pred <- self$predict(data,errors=T)
    # Randomly sample and inverse transform
    self$transform[[2]](rnorm(samples,pred$expect,pred$sd))
  }
  
  self
}
