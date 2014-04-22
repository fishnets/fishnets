# FAO 3 letter species codes
#
# Copyright 2003-2013 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# Soundtrack:
# Notes:
# TODO Best way of accessing info by code, name, etc.
# TODO Add fisbase code and method to open summary page

Fao3a <- object('Fao3a')

#' Create a data.frame of FAO3a codes and save it to file
#'
#' @param directory Directory where data will be saved
Fao3a$create <- function(directory="."){

	# GET zip file from FAO website
	temp <- tempfile()
	download.file("ftp://ftp.fao.org/FI/STAT/DATA/ASFIS_sp.zip", temp)

	# FIND txt version in zip
	dfiles <- unzip(temp, list=TRUE)
	file <- dfiles[grep(".txt", dfiles$Name), 'Name']

	# LOAD txt version of table
	fao <- read.table(unz(temp, file), sep="\t", header=TRUE)
	unlink(temp)

	# LOAD from local extracted copy
	fao <- read.table("ASFIS_sp_Feb_2013.txt", sep="\t", header=TRUE)

	# FIX Order
	Order <- as.character(fao$Order)
	NewOrder <- paste(toupper(substring(Order, 1, 1)),
		tolower(substring(Order, 2)), sep="")

	# CREATE data.frame
	# code3a, name: character
	# family, order: factor
	fao3a <- data.frame(code3a=as.character(fao$X3A_CODE),
		name=as.character(fao$Scientific_name), family=fao$Family,
		order=as.factor(NewOrder), stringsAsFactors=FALSE)

	# SAVE to file
	save(fb,file=file.path(directory,"fao3a.RData"))
}

#' Read data from disk
#' 
#' @param directory Directory where data will be read from
Fao3a$read <- function(directory="."){
  load(file.path(directory,"fao3a.RData"))  
  fao3a
}

#' Some tests
Fao3a$tests <- function(){
  fao3a <- Fao3a$read()

	attach(fao3a)

	name[code3a == 'SKJ']
	family[code3a == 'SKJ']
	order[code3a == 'SKJ']

	detach(fao3a)
}

