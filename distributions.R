Distribution <- function(rfunc,...){
  self <- object('Dist')
  
  self$random <- function(){
    rfunc(1,...)
  }
  
  self
}

Fixed <- function(value){
  self <- extend(Distribution,'Fixed')
  
  self$random <- function(){
    value
  }
  
  self
}
Uniform <- function(...) extend(Distribution,'Uniform',runif,...)
Normal <- function(...) extend(Distribution,'Normal',rnorm,...)
Lognormal <- function(...) extend(Distribution,'Lognormal',rlnorm,...)
Beta <- function(...) extend(Distribution,'Beta',rbeta,...)

dists <- function(...){
  dists <- lapply(as.list(match.call(expand.dots=T)[-1]),eval)
  class(dists) <- c('Distributions','list')
  dists
}
