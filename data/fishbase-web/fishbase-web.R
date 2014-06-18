#' R script for creating and reading FishbaseWeb data
require(plyr)
require(stringr)
require(rfishbase)
require(XML)
require(RCurl)

FishbaseWeb <- object('FishbaseWeb')

#' Get a table from a Fishbase URL
#' 
#' @param url URL for the page
#' @param which Which table (defaults to the first)
FishbaseWeb$get_table <- function(url,which=1){
  while(T){
    table <- try(readHTMLTable(addr))
    if(length(table) != 0) break
  }
  if(is(table, "try-error")) return(NULL)
  else table[[which]]
}

#' Get "List of Population Characteristics records"
#' 
#' @param id Species id
FishbaseWeb$get_popcharlist <- function(id){
  data <- FishbaseWeb$get_table(
    paste0("http://www.fishbase.org/PopDyn/PopCharList.php?ID=",id)
  )
}

#' A function to format the `res` list created by `getFB.R` and
#' transform it into a data.frame so it can be used for fitting
#' Fishnets
#' 
FishbaseWeb$format <- function(){

	# Load the downloade FishBase data
	load("FB.RData")

	# Some functions used below to convert factors to numerics and normalis sex codes etc
	char <- function(x) gsub("^\\s+|\\s+$", "", as.character(x))
	num <- function(x) as.numeric(char(x))
	sex <- function(x){
	  x = char(x)
	  x[x=='male'] = 'M'
	  x[x=='female'] = 'F'
	  x[x=='unsexed'] = 'U'
	  x[x=='mixed'] = 'U'
	  x[x=='juvenile'] = 'U'
	  x[is.na(x) | !(x %in% c('M','F','U'))] = 'U'
	  factor(x,levels=c('M','F','U'))
	}

	# For each species in the `res` list create a data.frame
	item <- 0
	dflist <- lapply(res,function(data){
	  
	  # Useful information for debugging
	  item <<- item + 1
	  cat(item,data$info$name,"\n")
	  
	  # $info - taxomomic and swim mode information
	  # A few species don't have valid $info so catch that situation...
	  result = tryCatch(with(data$info,data.frame(
	    species = name,
	    genus = genus,
	    family = family,
	    order = order,
	    swimmode = tolower(char(swimMode))
	  )),error=function(error)error)
	  if("error" %in% class(result)){
	    # Something wrong with basic information for this species (for some strange reason) so just return an empty data.frame
	    return(data.frame())
	  }
	  
	  # $growth
	  #   data.frame: Loo(cm), Length Type, K(1/y), to(years), Sex, M(1/y), Temp C, Lm, Ø', Country, Locality, Questionable, Captive, id
	  # The $growth table contains the 'primary' data for each species
	  # All other data for the species is replicated for each row in growth
	  
	  # First subset $growth so that it only includes useful rows
	  growth = subset(data$growth,Captive=="No" & Questionable=="No")
	  if(nrow(growth)>0){
	    # Put wanted values into a result data.frame
	    growth = data.frame(
	      sex = sex(growth[,'Sex']),
	      ltype = char(growth[,'Length Type']),
	      linf = num(growth[,'Loo(cm)']),
	      k = num(growth[,'K(1/y)']),
	      t0 = num(growth[,'to(years)']),
	      m = num(growth[,'M(1/y)']),
	      temp = num(growth[,7]),
	      country = char(growth[,'Country'])
	    )
      # If t0 is missing assume 0
	    growth$t0[is.na(growth$t0)] = 0
	    # Merge with info, keeping all growth estimates
	    result = merge(result,growth,all.y=T)
	  } else {
	    # No valid growth/,ortality daa for this species so return empty data.frame
	    return(data.frame())
	  }
	    
	  # $maturity
	  #   data.frame: Lm(cm), Length(cm)Low, Length(cm)Hig, Age(y)Low, Age(y)Hig, tm(y), Sex, Country, Locality, LengthType, id
	  maturity = data$maturity
	  if(is.data.frame(maturity)){
	    # Normalise sex
	    maturity$Sex = sex(maturity$Sex)
	    # Summarise by sex and country
	    maturity = ddply(maturity,.(Country,Sex),function(subset){
	      data.frame(
	        lmat = median(subset[,'Lm(cm)'],na.rm=T),
	        lmatmin = median(num(subset[,'Length(cm)Low']),na.rm=T),
	        lmatmax = median(num(subset[,'Length(cm)Hig']),na.rm=T),
	        amat = median(num(subset[,'tm(y)']),na.rm=T),
	        amatmin = median(num(subset[,'Age(y)Low']),na.rm=T),
	        amatmax = median(num(subset[,'Age(y)Hig']),na.rm=T)
	      )
	    })
	    names(maturity)[1:2] = c('country','sex')
	    result = merge(result,maturity,all.x=T,by=c('country','sex'))
	  }
	  
	  # $lw
	  #   data.frame:Score, a, b, Doubtful?, Sex, Length(cm), Length type, r2, SD, b, SD, log10, a, n, Country,Locality, id
	  lw = data$lw
	  if(is.data.frame(lw)){
	    # Use the "Doubtful?" column to restrict data
	    lw = lw[!grepl("Yes|yes",(lw[,4])),]
	    if(nrow(lw)>0){
	      # Normalise sex
	      lw$Sex = sex(lw$Sex)
	      # Summarise by sex and country
	      lw = ddply(lw,.(Country,Sex),function(subset){
	        data.frame(
	          a = median(num(subset[,'a']),na.rm=T),
	          b = median(num(subset[,'b']),na.rm=T)
	        )
	      })
	      names(lw)[1:2] = c('country','sex')
	      result = merge(result,lw,all.x=T,c('country','sex'))
	    }
	  }
	  
	  # $ecology
	  #   data.frame: tlevel, feedingHabit, feedingType, id
	  # This appear to have only one row, so merge that in
	  ecology = data$ecology
	  if(is.data.frame(ecology)){
	    if(nrow(ecology)>0){
	      result$trophic = num(ecology$tlevel[1])
	      result$diet = char(ecology$feedingType[1])
	      result$feeding = char(ecology$feedingHabit[1])
	      result$feeding[result$feeding==""] <- NA
	    }
	  }
	  
	  # $habitat
	  #   data.frame: habitat, id
	  # This data.frame can have as few as 2 rows and as many as 5 (or more?) depending on
	  # what habitat information is available.
	  habitat = data$habitat
	  if(is.data.frame(habitat)){
	    if(nrow(habitat)>0){
	      # First row is "habit" (e.g. benthopelagic)
	      value <- data$habitat[1,1]
	      result$habit <- if(is.null(value)) NA else char(value)
        # Second row is "migration" (e.g. anadromous) but can have additional information
        # on depth after the comma - discard the latter
        value <- data$habitat[2,1]
        result$migration <- if(is.null(value)) NA else strsplit(char(value),",")[[1]][1]
        # The "domain" information can get mixed in here so remove it
        result$migration[result$migation %in% c('marine','freshwater','brackish')] <- NA
        # Depth information can be on any of the pther rows. So paste all the rows together and then
        # grep it for depths
        value <- paste(data$habitat[,1],collapse=",")
        value <- str_match(value,"depth range (\\d+)*\\?* - (\\d+)*\\?* m")[1,2:3]
        result$depthmin <- num(value[1])
	      result$depthmax <- num(value[2])
	    }
	  }
	  
	  # $reproduction
	  #   data.frame: mode, fertilization, frequency, batch
	  # This appear to have only one row, so merge that in
    # This data appears to be incorrect see
    # So don't do this...
	  #reproduction = data$reproduction
	  #if(is.data.frame(reproduction)){
	  #   if(nrow(reproduction)>0){
	  #    result$repro = reproduction$mode[1]
	  #    result$fertil = reproduction$fertilization[1]
	  #    result$spawnfreq = reproduction$frequency[1]
	  #    result$spawnbatch = reproduction$batch[1]
	  #  }
	  #}
	  
	  # $fecundity
	  #   data.frame: country, location, AFmin, AFmax, RFmin, RFmean, RFmax, a, b
	  # AF = absolute fecundity
	  # RF = relative fecundity
	  # A very casula look at the data suggests that AF is more commonly available
	  # Here, rather than attempt to merge by country, just take means
	  fecundity = data$fecundity
	  if(is.data.frame(fecundity)){
	    if(nrow(fecundity)>0){
	      result$fecundmin = median(fecundity$AFmin,na.rm=T)
	      result$fecundmax = median(fecundity$AFmax,na.rm=T)
        # Sometimes AFMin>AFMax, sometimes both are 0, so create a new variable
	      result$fecundity = max(1,fecundity$AFmin,fecundity$AFmax,na.rm=T)
	    }
	  }
	  
	  result
	})

	# Row bind all the data.frames together
	fb <- do.call(rbind.fill,dflist)

	#' Augment the downloaded data with missing fields
	
	#' Currently, the downloaded data does not have class names. Currently, a lot of the other
	#' object in Fishnets assume this to be in the data that is fitted to (e.g TaxomicImputer)
	#' So we add class in here using an order to class lookup table generated from the Fishbase200 data
	#' (which does have class in it) using ClassLookupper:
	#' 		> cl <- ClassLookupper()
	#' 		> cl$fit(Fishbase2000$read('data/fishbase-2000'))
	#' 		> write.table(cl$table,file='data/fishbase-web/classes.txt',row.names=F,quote=F,sep='\t')
	#' This table was then augmented with order-class pairs that were missing
	fb <- merge(fb,read.table('classes.txt',header=T),all.x=T)

	#' Currently, the downloaded data does not have lmax and amax from the "http://www.fishbase.de/PopDyn/PopCharList.php?ID=" page
	#' So temporarily add the values for each species in Fishbase2000
	#		> fb <- Fishbase2000$read('data/fishbase-2000')
	#' 		> data <- ddply(fb,.(species),summarise,lmax=median(lmax,na.rm=T),amax=median(tmax,na.rm=T))
	#'		> write.table(data,file='data/fishbase-web/popchar.txt',row.names=F,quote=F,sep='\t')
	fb <- merge(fb,read.table('popchar.txt',header=T,sep='\t'),all.x=T)

  # Create factors
  fb <- within(fb,{
    swimmode <- factor(swimmode)
    diet <- factor(diet)
    feeding <- factor(feeding)
    habit <- factor(habit)
    migration <- factor(migration)
  })
  
	# Save
	save(fb,file="fishbase-web.RData")

}

#' Load from disk and return the FishbaseWeb data
FishbaseWeb$read <- function(directory='.'){
  load(file.path(directory,"fishbase-web.RData"))  
  fb
}

