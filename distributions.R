#' @name Distributions
#' 
#' Various probability distributions used to specify parameter priors
#' upon which to base sampling
#' 
#' @author Nokome Bentley

require(trapezoid)

Distribution <- function(rfunc,...){
  self <- object('Dist')

  self$random <- function(samples=1){
    rfunc(samples,...)
  }
  
  self
}

Fixed <- function(value){
  self <- extend(Distribution,'Fixed')
  
  self$random <- function(samples=1){
    rep(value,samples)
  }
  
  self
}
Uniform <- function(...) extend(Distribution,'Uniform',runif,...)
Normal <- function(...) extend(Distribution,'Normal',rnorm,...)
Lognormal <- function(...) extend(Distribution,'Lognormal',rlnorm,...)
Beta <- function(...) extend(Distribution,'Beta',rbeta,...)
Trapezoid <- function(...) extend(Distribution,'Trapezoid',rtrapezoid,...)
Triangle <- function(min,mode,max) extend(Distribution,'Triangle',rtrapezoid,min=min,mode1=mode,mode2=mode,max=max)

#' Create a list of distributions
#' The list is assigned a class so that a Fishnet knows how to deal with it
dists <- function(...){
  dists <- lapply(as.list(match.call(expand.dots=T)[-1]),eval)
  class(dists) <- c('Distributions','list')
  dists
}
