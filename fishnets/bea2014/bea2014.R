#' Development of a Fishnet for Bentley et al 2014
#' 
#' @author Nokome Bentley

# Source in the package
source('collate.R')

# Load the fishbase data
fb <- fishbase2000$read('data/fishbase-2000')

# Train some Fishnet nodes with the training set

tuning <- list()

variable = 'k'
tuning[[variable]] <- list()

# Taxonomic imputer node
ti <- Taxonomic.imputer(variable,c(log,exp),1)
tuning[[variable]][['ti']] <- ti$tune(fb)

# Glmer
glm <- Glmer(log(k)~class,exp,3,FALSE)
glm$tune(fb,c(
  #log(k)~class,
  #log(k)~class+order,
  #log(k)~class+order+family,
  #log(k)~class+order+family+genus,
  #log(k)~class+order+family+genus+species,
  #log(k)~class+order+family+genus+species+log(linf),
  log(k)~class+order+family+genus+species+log(linf)+bodyshape1
))

svm <- Svmer(log(k)~log(linf)+class+order+family+genus+species+foodtroph,exp)
svm$cross(fb,10)




# For the paper we examine out-of-sample predictive ability
# of the bea2014 Fishnet for three species:
#   Atlantic cod 'Gadus morhua'
#   Skipjack tuna 'Katsuwonus pelamis'
#   Bluenose 'Hyperoglyphe antarctica'

# Create input data.frames for each of the test species
cod <- data.frame(
  species = 'Gadus morhua',
  genus = 'Gadus',
  family = 'Gadidae',
  order = 'Gadiformes',
  class = 'Actinopterygii'
)
skj <- data.frame(
  species = 'Katsuwonus pelamis',
  genus = 'Katsuwonus',
  family = 'Scombridae',
  order = 'Perciformes',
  class = 'Actinopterygii'
)
bns <- data.frame(
  species = 'Hyperoglyphe antarctica',
  genus = 'Hyperoglyphe',
  family = 'Centrolophidae',
  order = 'Perciformes',
  class = 'Actinopterygii'
)

# Create a training data set which excludes those species
test_species = c('Gadus morhua','Katsuwonus pelamis','Hyperoglyphe antarctica')
fb_train <- subset(fb,!(species %in% test_species))

