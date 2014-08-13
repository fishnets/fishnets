
saver <- function(x, ..., name, path='C:/PROJECTS/FISHNETS/res/') {
  save(x, ..., file=paste(path,name,'.Rdata',sep=''))
}

pdfr <- function(x, ..., name, path='C:/PROJECTS/FISHNETS/res/') {
  pdf(..., file=paste(path,name,'.pdf',sep=''))
  print(x)
  dev.off()
}
