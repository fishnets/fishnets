#' fishnets uses an object oriented approach based on R environments
#' It is very similar to the approach described at http://www.lemnica.com/esotericR/Introducing-Closures/
#' and used in the proto package (http://cran.r-project.org/web/packages/proto/index.html)
#' Rather than introduce depenencies like proto, we just define this simple convieience function.
object <- function(class.name){
  self <- new.env()
  class(self) <- class.name
  self
}

extend <- function(base,class.name){
  self <- base()
  class(self) <- class.name
  self
}