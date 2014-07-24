#' A `Node` for matural mortality based on
#' [Charnov et al 2013]() Equation 3
MCharnovEtAl2013Fitted <- function(){
  self <- extend(MCharnovEtAl2013,'MCharnovEtAl2013Fitted')

  self$fit <- function(data){
    self$glm <- glm(log(m)~log(lmat/linf)+log(k),data)
  }
  
  self$predict <- function(data){
    predict.glm(self$glm,newdata=data,type='response',se.fit=TRUE)
  }
  
  self
}
