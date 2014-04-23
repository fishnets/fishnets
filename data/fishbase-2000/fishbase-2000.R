#' R script for processing FishBase data from CD-ROM circa 2000

Fishbase2000 <- object('Fishbase2000')

#' Create and write to disk a data.frame of Fishbase 2000 data
Fishbase2000$create <- function(directory){
  ## Read in the data
  for(name in c('ECOLOGY','FAMILIES','MORPHDAT','POPCHAR','POPGROWTH','POPLW','REPRODUC','SPAWNING','SPECIES','STOCKS')){
  	data = read.table(paste(name,'.txt',sep=''),header=T,colClasses="character",sep=",")
  	assign(name,data)
  }
  
  ##########################################################################################
  ## For each table, check, clean up and normalise data...
  
  ## ECOLOGY
  names(ECOLOGY)
  #[1] "StockCode"      "SpecCode"       "Stream"         "Lakes"          "Cave"          
  #[6] "Cave2"          "Estuaries"      "Mangroves"      "Intertidal"     "Soft"          
  #[11] "Rocky"          "Marine"         "Oceanic"        "Neritic"        "Coral.reefs"   
  #[16] "Soft.bottom"    "Hard.bottom"    "Sea.grass.beds" "Macrophyte"     "Main.food"     
  #[21] "Herbivory2"     "Feeding.Type"   "EcoTroph"       "DietTroph"      "FoodTroph" 
  unique(ECOLOGY$Main.food)
  
  ## FAMILIES
  names(FAMILIES) 
  # "FamCode"    "Family"     "CommonName" "Order"      "Class"     "ReprGuild"  "SwimMode"
  unique(FAMILIES$Family) #OK 
  unique(FAMILIES$Order) 
  # Get rid of some unecessary text...
  FAMILIES$Order[FAMILIES$Order=='Siluriformes (catfish)'] = 'Siluriformes'
  FAMILIES$Order[FAMILIES$Order=='Characiformes (characins)'] = 'Characiformes'
  # CommonName
  # Unecessary for classification but potenitally useful for interpretation so retain
  # Class
  unique(FAMILIES$Class)
  # Remove unecessary paranthetical text
  FAMILIES$Class = sub(' \\(.*\\)','',FAMILIES$Class)
  unique(FAMILIES$ReprGuild) 
  # Just make lower case
  FAMILIES$ReprGuild = tolower(FAMILIES$ReprGuild)
  unique(FAMILIES$SwimMode)
  
  ## SPECIES
  names(SPECIES) 
  #"SpecCode" "Genus"    "Species"  "FamCode" 
  unique(SPECIES$Genus)
  unique(SPECIES$Species)
  
  ## STOCKS
  names(STOCKS) 
  # "StockCode" "SpecCode"  "Level"     "Northern"  "Southern"  "TempMin"  "TempMax"
  unique(STOCKS$Level) #OK
  # The following are numeric...
  STOCKS$Northern = as.numeric(STOCKS$Northern)
  STOCKS$Southern = as.numeric(STOCKS$Southern)
  STOCKS$TempMin = as.numeric(STOCKS$TempMin)
  STOCKS$TempMax = as.numeric(STOCKS$TempMax)
  
  ## POPCHAR
  names(POPCHAR) # "Stockcode" "Speccode"  "Sex"       "Wmax"      "Lmax"      "tmax"
  names(POPCHAR) = c( "StockCode","SpecCode","Sex","Wmax","Lmax","tmax") #Make names consistent
  unique(POPCHAR$Sex)
  POPCHAR$Wmax = as.numeric(POPCHAR$Wmax)
  POPCHAR$Lmax = as.numeric(POPCHAR$Lmax)
  POPCHAR$tmax = as.numeric(POPCHAR$tmax)
  
  ## POPGROWTH
  names(POPGROWTH) 
  #[1] "StockCode"     "SpecCode"      "Sex"           "Loo"           "TLinfinity"    "K"             "to"            "Type"         
  #[9] "Winfinity"     "LinfLmax"      "Auxim"         "LogKLogLoo"    "tmax"          "M"             "MethodM"       "Mquality"     
  #[17] "Lm"            "LmLoo"         "LmSex"         "TypeLm"        "LmMale"        "LmLooMale"     "LmFemale"      "LmLooFemale"  
  #[25] "GrowthEnviron" "Temperature"   "DeltaT"
  
  unique(POPGROWTH$Sex)
  # Merge mixed and unsexed levels into one
  POPGROWTH$Sex[POPGROWTH$Sex=="mixed"] = "unsexed"
  
  POPGROWTH$K = as.numeric(POPGROWTH$K)
  POPGROWTH$to = as.numeric(POPGROWTH$to)
  unique(POPGROWTH$Type) #Fork length/ standard length/ total length etc
  POPGROWTH$Winfinity = as.numeric(POPGROWTH$Winfinity)
  POPGROWTH$LinfLmax = as.numeric(POPGROWTH$LinfLmax)
  unique(POPGROWTH$Auxim) #doubtful, unchecked, within ellipse
  POPGROWTH$LogKLogLoo = as.numeric(POPGROWTH$LogKLogLoo)
  POPGROWTH$tmax = as.numeric(POPGROWTH$tmax)
  POPGROWTH$M = as.numeric(POPGROWTH$M)
  unique(POPGROWTH$MethodM) #from tmax & Hoenig's model, tagging-recapture data etc
  unique(POPGROWTH$Mquality) #0 or 1
  unique(POPGROWTH$GrowthEnviron) # captivity,open waters
  POPGROWTH$Temperature = as.numeric(POPGROWTH$Temperature)
  POPGROWTH$DeltaT = as.numeric(POPGROWTH$DeltaT)
  
  ##POPLW
  names(POPLW) 
  # "StockCode"   "SpecCode"    "Type"        "LmaxCompare" "EsQ"        "Number"      "Sex"         "a"           "b" 
  unique(POPLW$Type)
  
  unique(POPLW$Sex)
  # Merge mixed and unsexed levels into one
  POPLW$Sex[POPLW$Sex=="mixed"] = "unsexed"
  
  unique(POPLW$Type)
  POPLW$a = as.numeric(POPLW$a)
  POPLW$b = as.numeric(POPLW$b)
  
  ##REPRODUC
  names(REPRODUC)
  #"StockCode"     "SpecCode"      "ReproMode"     "Fertilization" "Spawning"      "Batch.spawner" "RepGuild1"     "RepGuild2" 
  names(REPRODUC)[match("Batch.spawner",names(REPRODUC))] = c("BatchSpawner")
  unique(REPRODUC$ReproMode) # dioecism, parthenogenesis etc 
  unique(REPRODUC$Fertilization) # external, in mouth etc 
  unique(REPRODUC$Spawning)
  REPRODUC$Spawning = tolower(REPRODUC$Spawning) #  once in a lifetime, seasonal  etc
  unique(REPRODUC$BatchSpawner) # 0 or 1
  unique(REPRODUC$RepGuild1)
  REPRODUC$RepGuild1 = tolower(REPRODUC$RepGuild1) #  bearers,guarders etc
  unique(REPRODUC$RepGuild2)
  REPRODUC$RepGuild2 = tolower(REPRODUC$RepGuild2) #  nesters, brood hiders
  
  ##MORPHDAT
  names(MORPHDAT)
  #[1] "StockCode"           "Speccode"            "Body.Shape.I"        "Body.Shape.II"       "Forehead"           
  #[6] "OperculumPresent"    "Pos.of.mouth"        "Mandible.teeth"      "Mandible.teeth.T1"   "Mandible.teeth.T2"  
  #[11] "Mandible.rows.min"   "Mandible.rows.max"   "Maxilla.teeth"       "Maxilla.teeth.T1"    "Maxilla.teeth.T2"   
  #[16] "Maxilla.rows.min"    "Maxilla.rows.max"    "Vomerine.teeth"      "Vomerine.teeth.T1"   "Vomerine.teeth.T2"  
  #[21] "Vomerine.rows.min"   "Vomerine.rows.max"   "Palatine"            "Palatine.teeth.T1"   "Palatine.teeth.T2"  
  #[26] "Palatine.rows.min"   "Palatine.rows.max"   "Pharyngeal.teeth"    "Pharyngeal.teeth.T1" "Pharyngeal.teeth.T2"
  #[31] "Pharyngeal.rows.min" "Pharyngeal.rows.max" "Teeth.on.tongue"     "Lipsteeth"           "Dermal.teeth"       
  #[36] "Type.of.scales"      "Scutes"              "Keels"               "Maximum.depth" 
  # Give more consistent names..
  names(MORPHDAT)[1:4] = c('StockCode','SpecCode','BodyShape1','BodyShape2')
  # Just focus on the two body shape variables..
  MORPHDAT = MORPHDAT[,c('StockCode','SpecCode','BodyShape1','BodyShape2')]
  
  ##########################################################################################
  
  # Merge everything into a common dataframe
  # Starting with POPGROWTH
  fb = POPGROWTH
  nrow(fb)
  
  fb = merge(fb,
    SPECIES,
    by="SpecCode",all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    FAMILIES,
    by="FamCode",all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    STOCKS,
    by=c("SpecCode","StockCode"),all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    aggregate(POPCHAR[,c('Wmax','Lmax','tmax')],with(POPCHAR,list(SpecCode=SpecCode,StockCode=StockCode,Sex=Sex)),mean,na.rm=T),
    by=c("SpecCode","StockCode","Sex"),all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    aggregate(POPLW[,c('a','b')],with(POPLW,list(SpecCode=SpecCode,StockCode=StockCode,Sex=Sex)),mean,na.rm=T),
    by=c("SpecCode","StockCode","Sex"),all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    REPRODUC,
    by=c("SpecCode","StockCode"),all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    MORPHDAT,
    by=c("SpecCode","StockCode"),all.x=T
  )
  nrow(fb)
  
  fb = merge(fb,
    ECOLOGY,
    by=c("SpecCode","StockCode"),all.x=T
  )
  nrow(fb)
  
  ##########################################################################################
  # Tidy some things a little...
  
  # There are two tmax fields (one from POPCHAR, one from POPGROWTH). Just use one...
  sum(!is.na(fb$tmax.x)) #125
  sum(!is.na(fb$tmax.y)) #1924
  plot(tmax.y~tmax.x,fb); abline(a=0,b=1)
  fb$tmax = fb$tmax.x
  fb$tmax.x = NULL
  fb$tmax.y = NULL
  
  # Drop estimates from in captivity
  table(fb$GrowthEnviron) #228 versus 4973
  fb = subset(fb,GrowthEnviron!="captivity")
  fb$GrowthEnviron <- NULL #All just one value, so set to null
  
  ##Make NA any M estimates that are "from tmax & Hoenig's model"
  fb$M[fb$MethodM=="from tmax & Hoenig's model"] = NA
  
  # For many factors, there are blank, whitespace strings instead of NAs. 
  # make these blanks NAs
  for(name in names(fb)){
    if(is.character(fb[,name])){
      blanks = (!is.na(fb[,name])) & nchar(gsub("\\s","",fb[,name]))==0
      fb[blanks,name] = NA
    }
  }
  
  # Turn variables into factors, numeric etc, do some renaming and
  # remove variables that are not needed
  fb = within(fb,{
    SpecCode = as.factor(SpecCode)
    StockCode = as.factor(StockCode)
    Sex = as.factor(Sex)
    FamCode = as.factor(FamCode)
    Loo = as.numeric(Loo)
    LooTL = as.numeric(TLinfinity); TLinfinity <- NULL
    K = as.numeric(K)
    to = as.numeric(to)
    LType = as.factor(Type); Type <- NULL
    Woo = as.numeric(Winfinity); Winfinity <- NULL
    LinfLmax <- NULL
    Auxim <- NULL
    LogKLogLoo <- NULL
    Mmethod <- as.factor(MethodM); MethodM <- NULL
    Mquality <- as.factor(Mquality)
    Lm = as.numeric(Lm)
    LmLoo <- NULL
    LmSex <- NULL
    LmType <- as.factor(TypeLm); TypeLm <- NULL
    LmMale <- NULL
    LmLooMale <- NULL
    LmFemale <- NULL
    LmLooFemale <- NULL
    Temperature <- as.numeric(Temperature)
    TemperatureDelta <- as.numeric(DeltaT); DeltaT <- NULL
    Genus <- as.factor(Genus)
    Species <- as.factor(Species)
    Family <- as.factor(Family)
    Order <- as.factor(Order)
    Class <- as.factor(Class)
    ReprGuild <- as.factor(ReprGuild)
    SwimMode <- as.factor(SwimMode)
    Level <- as.factor(Level)
    
    ReproMode <- as.factor(ReproMode)
    Fertilization <- as.factor(Fertilization)
    Spawning <- as.factor(Spawning)
    BatchSpawner <- as.factor(BatchSpawner)
    RepGuild1 <- as.factor(RepGuild1)
    RepGuild2 <- as.factor(RepGuild2)
    BodyShape1 <- as.factor(BodyShape1)
    BodyShape2 <- as.factor(BodyShape2)
  })
  for(name in c(
    'Stream','Lakes','Cave','Cave2','Estuaries','Mangroves','Intertidal',
    'Soft','Rocky','Marine','Oceanic','Neritic','Coral.reefs','Soft.bottom','Hard.bottom',
    'Sea.grass.beds','Macrophyte','Main.food','Herbivory2','Feeding.Type'
  )) fb[,name] = as.factor(fb[,name])
  for(name in c('EcoTroph','DietTroph','FoodTroph')) fb[,name] = as.numeric(fb[,name])
  
  # Make all names lower case and change some to reduce ambiguities
  names(fb) = tolower(names(fb))
  names(fb)[names(fb)=='lm'] = 'lmat'
  names(fb)[names(fb)=='loo'] = 'linf'
  
  # Make species a combination of genus and species names
  fb$species = factor(with(fb,paste(genus,species)))
  
  save(fb,file=file.path(directory,"fishbase-2000.RData"))
}

#' Load from disk and return the Fishbase 2000 data
Fishbase2000$read <- function(directory){
  load(file.path(directory,"fishbase-2000.RData"))  
  fb
}
