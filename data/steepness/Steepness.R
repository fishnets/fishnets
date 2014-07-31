#' A `Node` for steepness based on
#' [Myers et al 1999](http://www.nrcresearchpress.com/doi/abs/10.1139/f99-201). 
#' Note that Myers et al. fit a Ricker model to the stock recruit data, but then calculate the Beverton-Holt steepness from the Ricker \alpha. This has been shown to be overly conservative, and the node will post a warning to this effect.

source('collate.R')

Steepness <- object('SteepnessMyersEtAl1999')

#' Read data from disk
#' 
#' @param directory Directory where data will be read from
Steepness$create <- function(directory = "."){
  
  steep <- read.table(file.path(directory,'Myers_et_al_steepness_extended.csv'),sep=',',comment.char = c(''), stringsAsFactors = F,header=T)
  steep
}
  
#' Merge steepness data with life-history data to inform an empirical node
#' 
#' @param steep Steepness predictions for a subset of spescies
#' @param lh_data Life-history data to append
#' @param merge_by column in life history data to merge by; defaults to 'species'

Steepness$merge <- function(steep,lh_data,merge_by='species'){
  
  # matches <- match(as.character(fb$species),steep$species)
  # lh_matches <- lh_data[!is.na(matches),]
  
  steep_merged <- merge(steep,lh_data,by.x='species',by.y=merge_by,all.y=T)
  steep_merged
  
}




