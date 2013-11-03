require(trapezoid)

Distribution <- function(rfunc,...){
  self <- object('Dist')
  
  self$random <- function(n=1){
    rfunc(n,...)
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
Trapezoid <- function(...) extend(Distribution,'Trapezoid',rtrapezoid,...)
Triangle <- function(min,mode,max) extend(Distribution,'Triangle',rtrapezoid,min=min,mode1=mode,mode2=mode,max=max)

dists <- function(...){
  dists <- lapply(as.list(match.call(expand.dots=T)[-1]),eval)
  class(dists) <- c('Distributions','list')
  dists
}
