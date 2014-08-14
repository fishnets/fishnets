
saver <- function(..., name, path='C:/PROJECTS/FISHNETS/res/') {
  save(..., file=paste(path,name,'.Rdata',sep=''))
}

loader <- function(name, path='C:/PROJECTS/FISHNETS/res/') {
  do.call(load,list(file=paste(path,name,'.Rdata',sep=''),envir=globalenv()))
}

pdfr <- function(x, ..., name, path='C:/PROJECTS/FISHNETS/res/') {
  pdf(..., file=paste(path,name,'.pdf',sep=''))
  print(x)
  dev.off()
}
