# bayesianMPRA
An analytical framework for MPRA data with two main goals:
1. Provide a fully Bayesian, generative statistical model for the raw count data
2. Integrate genomic annotations to provide an informative empirical prior to increase the statistical power of the analysis in an interpretable way

This package is still evolving and may be broken at times. If you run into a problem please create an issue on this github page.

## Installation  

### Dependencies  

The package requires the tidyverse packages, `rstan`, `coda`, and `fitdistrplus` and suggests the bioconductor package `preprocessCore`.

```
install.packages(pkgs = c('tidyverse', 'rstan', 'coda', 'fitdistrplus'))

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

#### Example usage

```
> mpra_data
# A tibble: 25,480 x 10
          snp_id allele K562_minP_DNA1 K562_minP_DNA2 K562_CTRL_minP_RNA1 K562_CTRL_minP_RNA2 K562_CTRL_minP_RNA3
           <chr>  <chr>          <int>          <int>               <int>               <int>               <int>
 1 1 3691528 2/3    ref           1318            885                1243                 593                1109
 2 1 3691528 2/3    ref           1133            731                 613                 395                1562
 3 1 3691528 2/3    ref            423            276                 192                  83                 257
 4 1 3691528 2/3    ref            459            324                 189                 108                 303
 5 1 3691528 2/3    ref            396            254                 325                  95                 324
 6 1 3691528 2/3    ref            905            555                 443                 255                1390
 7 1 3691528 2/3    ref            585            380                 363                 219                 813
 8 1 3691528 2/3    ref            525            385                 291                 124                 216
 9 1 3691528 2/3    ref            831            546                 479                 302                1472
10 1 3691528 2/3    ref            169            117                  55                  43                 281
# ... with 25,470 more rows, and 3 more variables: K562_CTRL_minP_RNA4 <int>, K562_CTRL_minP_RNA5 <int>,
#   K562_CTRL_minP_RNA6 <int>

> predictors
# A tibble: 910 x 3
             snp_id DeepSeaDnaase gkmerDelta
             <fctr>         <dbl>      <dbl>
 1 10 101274365 1/3     -0.131800   0.422375
 2 10 101275206 1/2     -0.295580  -0.664241
 3 10 101292390 2/3     -0.107650   2.206986
 4  10 45949254 2/3      0.090070  -1.792320
 5  10 45951081 1/2      0.263740  -0.970451
 6  10 45952929 1/2      0.206140   1.433608
 7  10 45958275 1/3     -0.406990  -3.264995
 8  10 45959385 2/3     -0.094865   2.011715
 9  10 45963777 2/3      0.554320   2.323613
10  10 45965085 1/3     -0.277360  -1.377249
# ... with 900 more rows

bayesian_mpra_analyze(mpra_data, predictorss, out_dir = '~/bayesianMPRA/analysis_outputs/', num_cores = 8)

```
#### Outputs  
The output will be a data frame with a row for each `snp_id` containing a nested column of the count data, along with columns giving the MLE negative binomial parameters, the estimated gamma prior, and a binary call of "functional" or "non-functional" based on the snp_id's posterior. For functional variants (according to a 95% HPD interval on transcriptional shift), the MCMC results, transcriptional shift HDI's and mean transcriptional shifts are written to .RData objects in `out_dir` with the appropriate snp_id.

## Notes

The `analysis` directory in this github repo is for the exploratory scripts used by the author to create this package. They are not installed with the package and should be ignored.
