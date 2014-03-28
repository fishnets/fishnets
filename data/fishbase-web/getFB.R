# getFB.R - DESC
# getFB.R

# Copyright 2013 Iago Mosqueira, JRC. Distributed under the GPL 2 or late.
# Maintainer: Iago Mosqueira, JRC
# Soundtrack:
# Notes:
# TODO: Add error checking on existence of tables at every step

library(rfishbase)
library(XML)
library(RCurl)

# UPDATE summaries
# updateCache()

# GET local fishbase
data(fishbase)

names(fish.data) <- unlist(lapply(fish.data, function(x) x$id))

# ids from web page http://www.fishbase.org/Topic/List.php?group=9
# vim :%s/.*\(ID=[0-9]*\).*/\1
# vim :v/ID=[0-9]*/d
# vim :%s/ID=//
ids <- sort(scan('ids.dat'))

# TEST
# ids <- ids[1:40]
# ids <- sample(ids, 15)

# CREATE output list
out <- vector(mode='list', length=length(ids))
names(out) <- ids

# LOOP {{{
for(i in ids){

	# RESTART
	if(!is.null(out[[as.character(i)]]))
		 next

	cat(i)

	# SPP
	spp <- vector(mode='list', length=8)
	names(spp) <- c('info', 'growth', 'maturity', 'lw', 'ecology', 'habitat',
		'reproduction', 'fecundity')

	# INFO {{{
	Info <- fish.data[[as.character(i)]]
  
  # TODO updatecache to get all names
  if(is.null(Info))
    Info <- list(id=i, Genus=NULL, Fanily=NULL, Order=NULL, ScientificName="NA NA")
  
	sppInfo <- list(id=i, genus=Info$Genus, family=Info$Family,
		order=Info$Order, name=Info$ScientificName, Species=strsplit(Info$ScientificName,
	 	" ")[[1]][2], swimMode="NA")

	# SWIM MODE
	addr <- paste0("http://www.fishbase.org/physiology/swimmodesummary.php?speccode=", i)
	tab <- try(readHTMLTable(addr))
  while(length(tab) == 0)
    tab <- try(readHTMLTable(addr))

	if(is(tab, "try-error")) {
		swimmode <- NULL
	} else {

		# TODO CHECK results, assumes space at start
		swimmode <- sub("^.", "", as.character(tab[[3]][1,4]))

		# No info
		if(swimmode == "")
			swimmode  <- "NA"
	}
	sppInfo$swimMode <- swimmode

	spp$info <- sppInfo # }}}
	
	# GROWTH {{{
	addr <- paste("http://www.fishbase.org/PopDyn/PopGrowthList.php?ID=",i, sep="")
	
	tab <- try(readHTMLTable(addr))
	
	if(is(tab, "try-error")) {
		growth <- NULL
	} else {

		# No growth info
		if(length(tab) == 0) {
			growth <- "NA"
			# BREAK if no growth info
			cat("[", i, "]")
			next
		} else {

			# DROP 1st column
			growth <- tab$dataTable[,-1]

			# ADD id
			growth <- cbind(growth, data.frame(id=i))
		}
	}
	spp$growth <- growth
	# }}}

	# MATURITY {{{
	addr <- paste("http://www.fishbase.org/Reproduction/MaturityList.php?ID=",i, sep="")

	tab <- try(readHTMLTable(addr))
	while(length(tab) == 0)
	  tab <- try(readHTMLTable(addr))
	
	if(is(tab, "try-error")) {
		mature <- NULL
	} else {
	
		if(!is.null(tab$dataTable)){
	
			# ADD names
			names(tab$dataTable) <- c('', 'Lm(cm)', 'Length(cm)Low','',
			'Length(cm)Hig', 'Age(y)Low','', 'Age(y)Hig','tm(y)','Sex','Country','Locality')

			# DROP empty columns
			mature <- tab$dataTable[,-c(1, 4, 7)]

			# SPLIT Lm(cm) column
			Lm <- strsplit(as.character(mature[,'Lm(cm)']), ' ')
			mature[,'Lm(cm)'] <- unlist(lapply(Lm, function(x) as.double(x[1])))
			mature$LengthType <- unlist(lapply(Lm, function(x) as.character(x[2])))
		
			# ADD id
			mature <- cbind(mature, data.frame(id=i))

		} else {
			mature <- "NA"
		}
	}
	spp$maturity <- mature
	# }}}

	# LENGTH-WEIGHT {{{
	addr <- paste("http://www.fishbase.org/PopDyn/LWRelationshipList.php?ID=",i, sep="")
	tab <- try(readHTMLTable(addr))
	while(length(tab) == 0)
	  tab <- try(readHTMLTable(addr))

	if(is(tab, "try-error")) {
		lw <- NULL
	} else if(is.null(tab[[3]])){
		lw <- NULL
	} else {
		names(tab)[3] <- "dataTable"
		
		lw <- tab$dataTable

		# FIX names
		names(lw)[c(6,7)] <- c("Length(cm)", "Length type")

		# ADD id
		lw <- cbind(lw, data.frame(id=i))

	}
	spp$lw <- lw
	# }}}

	# ECOLOGY {{{
	addr <- paste0("http://www.fishbase.org/summary/", i)

	page <- try(RCurl::getURLContent(addr))

	if(is(page, "try-error")) {
		ecology <- NULL
	} else {
		page <- strsplit(page, "\r\n")[[1]]

		# GET StockCode
		idx <- grep("FishEcologySummary", page, fixed=TRUE, value=TRUE)
		if(length(idx) == 0) {
			ecology <- NULL
		} else {

			mat <- regexpr("StockCode=[[:digit:]]*", idx)

		stkcd <- regmatches(idx, mat)[1]

		addr <- paste("http://www.fishbase.org/Ecology/FishEcologySummary.php?",
			stkcd, "&GenusName=", strsplit(Info$ScientificName, " ")[[1]][1],
			"&SpeciesName=", strsplit(Info$ScientificName, " ")[[1]][2], sep="")
    

		tab <- try(readHTMLTable(addr))
    while(length(tab) == 0)
      tab <- try(readHTMLTable(addr))

		# GET Feeding habit
		feedingHabit <- as.character(tab[[5]][2,2])

		# GET trophic level
		mat <- regexpr("StockCode=[[:digit:]]*", idx)

		if(!is(tab, "try-error"))
				if(length(tab) > 0)
					tlev <- suppressWarnings(as.double(as.character(tab[[6]][,2])))[4]

		# TODO: ADD feeding type
		feedingType <- ifelse(tlev < 2.20, 'primary',
			ifelse(tlev <= 2.8, 'omnivore', 'secondary'))
	
		# ADD id
		ecology <- data.frame(tlevel=tlev, feedingType=feedingType,
			feedingHabit=feedingHabit, id=i)
		}
	}
	spp$ecology <- ecology
	# }}}

	# HABITAT {{{

	# GET habitat from Info
	habitat <- unlist(strsplit(gsub(".Ref. [0-9]+.", "", Info$habitat), "; "))
	
	# ADD id
	if (!is.null(habitat))
		habitat <- data.frame(habitat=habitat, id=i)

	spp$habitat <- habitat
	# }}}
	
	# MORPH

	# POPCHAR: wmax lmax tmax

	# REPROD {{{
	addr <- paste0("http://www.fishbase.org/Reproduction/FishReproSummary.php?ID=", i, "&GenusName=", sppInfo$Genus, "&SpeciesName=", sppInfo$Species, "&StockCode=1")
	
	tab <- try(readHTMLTable(addr))
	while(length(tab) == 0)
	  tab <- try(readHTMLTable(addr))

	if(is(tab, "try-error")) {
		reproduction <- NULL
	} else {
		reproduction <- data.frame(
			mode=as.character(tab[[1]][1,2]),
			fertilization=as.character(tab[[1]][2,2]),
			frequency=as.character(tab[[1]][3,2]),
			batch=ifelse(grep("Yes", as.character(tab[[1]][4,2])) == 1, TRUE, FALSE))
	}
	spp$reproduction <- reproduction
	# }}}

	# FECUNDITY {{{
	addr <- paste0("http://www.fishbase.org/Reproduction/FecundityList.php?ID=", i, "&GenusName=", sppInfo$Genus, "&SpeciesName=", sppInfo$Species, "&StockCode=1")
	
	tab <- try(readHTMLTable(addr))
	while(length(tab) == 0)
	  tab <- try(readHTMLTable(addr))
	
	if(is(tab, "try-error")) {
		fecundity <- NULL
	} else {

		fecundity <- data.frame(country=as.character(tab[[2]][-1, 1]),
			location=as.character(tab[[2]][-1, 2]),
			AFmin=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 3]))),
			AFmax=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 4]))),
			RFmin=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 5]))),
			RFmean=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 6]))),
			RFmax=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 7]))),
			a=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 8]))),
			b=as.numeric(gsub(",", "", as.character(tab[[2]][-1, 9]))))
	}

	spp$fecundity <- fecundity
	# }}}

	# STORE spp
	out[[as.character(i)]] <- spp
  
	cat("\n")

} # }}}

# CHECKS for no growth data
idx <- unlist(lapply(out, function(x) is.null(x$growth)))
res <- out[!idx]

# ADDS FAO3A codes

# PERFECT matches
res <- lapply(res, function(x) {
	print(x$info$id)
	code <- fao3a[match(x$info$name, fao3a$name), 'code3a']
	if(is.na(code)) {
		# TRY dropping subspecies
		code <- fao3a[match(paste(strsplit(x$info$name, " ")[[1]][1:2], collapse=" "),
		fao3a$name), 'code3a']
		if(is.na(code)) {
			# TRY genus + spp
			code <- fao3a[match(paste(strsplit(x$info$name, " ")[[1]][1], "spp",
				collapse=" "), fao3a$name), 'code3a']
			if(is.na(code)) {
				code <- 'NA'
			}
		}
	}
	else
		x$info$fao3a <- code
	return(x)
})


# SAVE
save(out, res, file="FB.RData")
