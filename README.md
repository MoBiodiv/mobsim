# mobsim

## Spatial analysis of scale-dependent biodiversity changes ##

The package includes functions to simulated species distributions in space with controlled abundance distributions as well as controlled intraspecific aggregation. For analysis there are functions for species rarefaction and accumulation curves, species-area relationships,  endemics-area relationships and th distance-decay of community similarity.


## Installing mobsim

Note mobsim includes efficient C++ code. Therefore the package Rcpp and a
C++ compiler (e.g. Rtools for Windows) need to be installed first

[mobsim on GitHub](https://github.com/MoBiodiv/mobsim)

library(devtools)

install_github("MoBiodiv/mobsim", build_vignettes = T)

library(mobsim)

## Getting help

?mobsim

browseVignettes("mobsim")
