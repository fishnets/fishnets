## Fishnets : A library of multivariate priors for fish population dynamics parameters

### Status

This R package is currently in active development. Once the 
code interface stabilizes we hope to publish it on CRAN.

### Nodes

A `Fishnet` is made up of one or more nodes.
Each node is a predictor for a particular population dynamics parameter.
Fishnets has several generic `Node` classes i.e. those that can be applied to most, if not all, population dynamics parameters e.g. `Svmer`, `Glmer`. 
Fishnets also has some specific `Node` classes, i.e. those that apply to a specific population dynamics parameter e.g. `RecautoThorsonEtAl2014`, `RecsteepHeEtAl2006`.
See [issues](https://github.com/fishnets/fishnets/issues) for suggestions for more types of generic and specific nodes.

### Acknowledgements

This project has been supported by the [European Commision's Joint Research Centre](https://fishreg.jrc.ec.europa.eu/) and [Trident Systems](http://tridentsystems.co.nz/)
