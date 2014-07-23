#' A `Node` for matural mortality based on
#' [Charnov et al 2013]() Equation 3
MCharnovEtAl2013Fitted <- function(){
  self <- extend(MCharnovEtAl2013,'MCharnovEtAl2013Fitted')

  self$fit <- function(data){
    predicted <- log(self$predict(data))
    observed <- log(data$m)
    self$error <- sd(predicted-observed,na.rm=T)
  }
  
  self$predict <- function(data){
    with(data,((lmat/linf)^-1.5)*k)
  }
  
  self
}
