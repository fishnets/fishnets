#' Create an object
#' 
#' Fishnets uses an object oriented approach based on R environments.
#' It is very similar to the approach described at http://www.lemnica.com/esotericR/Introducing-Closures/
#' and used in the proto package (http://cran.r-project.org/web/packages/proto/index.html).
#' This function is a convienience function for creating a new environment ad setting its class
#'
#' @author Nokome Bentley
#' 
#' @param base Base object
#' @param class.name Class name to use
#' @param ... Other arguments to pass to base
object <- function(class.name){
  self <- new.env()
  class(self) <- class.name
  self
}

#' Extend an object (e.g. by adding new methods)
#' 
#' @author Nokome Bentley
#' 
#' @param base Base object
#' @param class.name Class name to use
#' @param ... Other arguments to pass to base
extend <- function(base,class.name,...){
  self <- base(...)
  class(self) <- c(class.name,class(self))
  self
}
