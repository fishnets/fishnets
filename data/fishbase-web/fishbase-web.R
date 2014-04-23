#' R script for creatring and reading in FishbaseWeb data
require(plyr)

FishbaseWeb <- object('FishbaseWeb')

#' A function to format the `res` list created by `getFB.R` and
#' transform it into a data.frame so it can be used for fitting
#' Fishnets
#' 
#' The does a minimum amound of normalisation, only what is necesary
#' to merge the data into a data.frame. Anoter script does normalisation to
#' prepare the data to be ready for analysis.
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
	    swimmode = swimMode
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
	        lmatlo = median(num(subset[,'Length(cm)Low']),na.rm=T),
	        lmathi = median(num(subset[,'Length(cm)Hig']),na.rm=T),
	        amat = median(num(subset[,'tm(y)']),na.rm=T),
	        amatlo = median(num(subset[,'Age(y)Low']),na.rm=T),
	        amathi = median(num(subset[,'Age(y)Hig']),na.rm=T)
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
	      result$trophic = ecology$tlevel[1]
	      result$feedtype = ecology$feedingType[1]
	      result$feedhabit = ecology$feedingHabit[1]
	    }
	  }
	  
	  # $habitat
	  #   data.frame: habitat, id
	  # This data.frame can have as few as 2 rows and as many as 5 (or more?) depending on
	  # what habitat information is available. Sometimes depth information is in the last rowm sometimes not.
	  # So, the approach taken here is to paste it all together and then grep it later
	  habitat = data$habitat
	  if(is.data.frame(habitat)){
	    if(nrow(habitat)>0){
	      result$habitat = paste(habitat[,1],collapse=",",sep="")
	    }
	  }
	  
	  # $reproduction
	  #   data.frame: mode, fertilization, frequency, batch
	  # This appear to have only one row, so merge that in
	  reproduction = data$reproduction
	  if(is.data.frame(reproduction)){
	    if(nrow(reproduction)>0){
	      result$repromode = reproduction$dioecism[1]
	      result$reprofertil = reproduction$fertilization[1]
	      result$reprofreq = reproduction$frequency[1]
	      result$reprobatch = reproduction$batch[1]
	    }
	  }
	  
	  # $fecundity
	  #   data.frame: country, location, AFmin, AFmax, RFmin, RFmean, RFmax, a, b
	  # AF = absolute fecundity
	  # RF = relative fecundity
	  # A very casula look at the data suggests that AF is more commonly available
	  # Here, rather than attempt to merge by country, just take means
	  fecundity = data$fecundity
	  if(is.data.frame(fecundity)){
	    if(nrow(fecundity)>0){
	      result$fecunditymin = median(fecundity$AFmin,na.rm=T)
	      result$fecunditymax = median(fecundity$AFmax,na.rm=T)
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

	# Save
	save(fb,file="fishbase-web.RData")

}

#' Load from disk and return the FishbaseWeb data
FishbaseWeb$read <- function(directory='.'){
  load(file.path(directory,"fishbase-web.RData"))  
  fb
}
