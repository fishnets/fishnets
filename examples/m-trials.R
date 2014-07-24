require(ggplot2)

# Currently, this script must be run in the Fishnets top level directory
# Source in the package
source('collate.R')

# Load the Fishbase data
fb <- FishbaseWeb$read('data/fishbase-web')

# Source in m relayted nodes (they might nt be in collate.R yet)
source('m-charnov-et-al-2013.R')
source('m-charnov-et-al-2013-fitted.R')

cea13 <- MCharnovEtAl2013()
cea13fit <- MCharnovEtAl2013Fitted()

cea13$cross(fb)

cea13fit
