#' A network node that predicts based on taxonomic group.
#' It moves up the taxonomic levels until it finds a level which has at least nmin
#' values for the variable and then uses the median (for numeric variables) or mode (for factors etc)
#' at that taxonomic level.
#' 
#' @author Nokome Bentley
#' 
#' @param predictand The name of the variable of being predicted
#' @param transform A pair of from/to transformation functions e.g. c(log,exp)
#' @param nmin Minimum sample size required for median/mode
TaxonomicImputer <- function(predictand,transform=c(identity,identity),nmin=1){
  self <- extend(Node,'TaxonomicImputer')
  
  self$predictand <- predictand
  self$transform <- transform
  self$predictors <- c('species','genus','family','order','class')
  self$nmin <- nmin
  
  self$fit <- function(data, ...){
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
          if(!is.null(defaults$n)){
            if(defaults$n>=self$nmin){
              best <- defaults
              break
            }
          }
        }
      }
      if(!is.null(best)){
        if(errors) return(best)
        else preds = c(preds,self$transform[[2]](best$expect))
      } else {
        if(errors) return(NA)
        else preds = c(preds,NA)
      }
    }
    preds
  }
  
  #' Tune nmin through simple search over a range of values
  self$tune <-function(data,trials=1:30){
    results <- NULL
    for(nmin in trials){
      cat("nmin:",nmin,"\n")
      self$nmin = nmin
      results <- rbind(results,data.frame(
        nmin = nmin,
        self$cross(data)$summary
      ))
    }
    list(
      best = results$nmin(which.min(results$mpe)),
      trials = results
    )
  }
  
  self$sample <- function(data){
    # Currently, only provides sampling for numeric predictands,
    # for factors just return the prediction
    if(!is.numeric(self$type())) return(self$predict(data))
    # Get prediction
    pred <- self$predict(data,errors=T)
    # Randomly sample and inverse transform
    self$transform[[2]](rnorm(nrow(data),mean=pred$expect,sd=pred$sd))
  }
  
  self
}
