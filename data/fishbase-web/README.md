% Extraction of data from fishbase
%
%

# TODO

- rfishbase::update.cache()
- Get final data.frame with all species
- Merge fao3a


# $info

- From pkg:::rfishbase
- character: id, genus, family, order, name

## swimMode

- $info['swimMode'], character

# $growth

- From http://www.fishbase.org/PopDyn/PopGrowthList.php?ID=id
- data.frame: Loo(cm), Length Type, K(1/y), to(years), Sex, M(1/y), Temp C, Lm, Ø', Country, Locality, Questionable, Captive

# $maturity
- From: http://www.fishbase.org/Reproduction/MaturityList.php?ID=id
- data.frame: Lm(cm), Length(cm)Low, Length(cm)Hig, Age(y)Low, Age(y)Hig, tm(y), Sex, Country, Locality, LengthType, id

# $lw

- From http://www.fishbase.org/PopDyn/LWRelationshipList.php?ID=id
- data.frame:Score, a, b, Doubtful?, Sex, Length(cm), Length type, r2, SD, b, SD, log10, a, n, Country,Locality, id

# $ecology

- From: http://www.fishbase.org/Ecology/FishEcologySummary.php?
- data.frame: tlevel, feedingHabit, feedingType, id
- feedingType
	- primary: 2.0 - 2.19
	- secondary: >=  2.8
	- omnivores: 2.2 - 2.79. 

# $habitat

- From pkg:::rfishbase
- data.frame: habitat, id

# $reproduction

- From: http://www.fishbase.de/Reproduction/FishReproSummary.php?ID=
- data.frame: mode, fertilization, frequency, batch

# $fecundity

- From: http://www.fishbase.de/Reproduction/FecundityList.php?ID=
- data.frame: country, location, AFmin, AFmax, RFmin, RFmean, RFmax, a, b

# FB data.frame

> names(fb)
 [1] "SpecCode"         "StockCode"        "Sex"              "FamCode"         
 [5] "Loo"              "K"                "to"               "M"               
 [9] "Mquality"         "Lm"               "Temperature"      "Genus"           
[13] "Species"          "Family"           "CommonName"       "Order"           
[17] "Class"            "ReprGuild"        "SwimMode"         "Level"           
[21] "Northern"         "Southern"         "TempMin"          "TempMax"         
[25] "Wmax"             "Lmax"             "a"                "b"               
[29] "ReproMode"        "Fertilization"    "Spawning"         "BatchSpawner"    
[33] "RepGuild1"        "RepGuild2"        "BodyShape1"       "BodyShape2"      
[37] "Stream"           "Lakes"            "Cave"             "Cave2"           
[41] "Estuaries"        "Mangroves"        "Intertidal"       "Soft"            
[45] "Rocky"            "Marine"           "Oceanic"          "Neritic"         
[49] "Coral.reefs"      "Soft.bottom"      "Hard.bottom"      "Sea.grass.beds"  
[53] "Macrophyte"       "Main.food"        "Herbivory2"       "Feeding.Type"    
[57] "EcoTroph"         "DietTroph"        "FoodTroph"        "tmax"            
[61] "TemperatureDelta" "LmType"           "Mmethod"          "Woo"             
[65] "LType"            "LooTL"

## df structure

- Keep all growth studies
- Lmax, tmax and wmax: one per spp/stock/sex (max)
- popchart: one per spp/stock/sex (mean)

## Convert to single df

http://r.789695.n4.nabble.com/Convert-quot-ragged-quot-list-to-matrix-td895283.html
