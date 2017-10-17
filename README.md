# bayesianMPRA
An analytical framework for MPRA data with two main goals:
1. Provide a fully Bayesian, generative statistical model for the raw count data with an empirical prior
2. Integrate genomic annotations to provide an informative empirical prior

This package is still evolving and may be broken at times. If you run into a problem please create an issue on this github page.

## Installation  

### Dependencies  

The package requires the tidyverse packages, `rstan` and suggests the bioconductor package `preprocessCore`.

```
install.packages(pkgs = c('tidyverse', 'rstan'))

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("preprocessCore")
```

### Package installation

Install the R package `devtools` if you don't have it already. Once you have that, the package is easily installed with:

```
devtools::install_github('andrewGhazi/bayesianMPRA') 
```
## Usage

Give the function `bayesian_mpra_analyze` a data frame of MPRA counts `mpra_data` and a data frame of functional `predictors` to perform the analysis. The outputs of functional variants are writen to `out_dir`. By default nonfunctional results are not written out because they can consume a lot of storage space. 

## Notes

The `analysis` directory in this github repo is for the exploratory scripts used by the author to create this package. They are not installed with the package and should be ignored.
