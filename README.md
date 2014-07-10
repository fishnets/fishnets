## Fishnets : A library of multivariate priors for fish population dynamics parameters

### Background

Fishnets is an R package for developing and using Bayesian networks for fish population dynamics parameters. 
See [Bentley (2014)](http://icesjms.oxfordjournals.org/content/early/2014/03/04/icesjms.fsu023.abstract) for a rationale and potential uses.

### Status

Fishnets is currently in active development. Once the code interface stabilizes we hope to publish it on CRAN.

### Nodes

A `Fishnet` is made up of one or more `Node`s.
Each node is a predictor for a particular population dynamics parameter.
Fishnets has several generic `Node` classes i.e. those that can be applied to most, if not all, population dynamics parameters e.g. `Svmer`, `Glmer`. 
Fishnets also has some specific `Node` classes, i.e. those that apply to a specific population dynamics parameter e.g. `RecautoThorsonEtAl2014`, `RecsteepHeEtAl2006`.
See [issues](https://github.com/fishnets/fishnets/issues) for suggestions for more types of generic and specific nodes.

### Acknowledgements

This project has been supported by the [European Commision's Joint Research Centre](https://fishreg.jrc.ec.europa.eu/) and [Trident Systems](http://tridentsystems.co.nz/). Many thanks to the editors and contributors of [FishBase](www.fishbase.org)(Froese & Pauly 2011) for providing a valuable source of information on population dynamics parameters

### References

Bentley, N. (2014). Data and time poverty in fisheries estimation: potential approaches and solutions. ICES Journal of Marine Science: Journal du Conseil, doi:10.1093/icesjms/fsu023.

Froese, R. and D. Pauly. Editors. (2011). FishBase. World Wide Web electronic publication. www.fishbase.org
