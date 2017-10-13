# let's lay out the skeleton for the package mostly adapted over from
# stanModelInformativePrior and changed to be generalized functions where
# applicable

# library(tidyverse)
# library(magrittr)
# library(rstan)
# library(parallel)
# library(fitdistrplus)
# source('R/mledistModified.R')

# ### Stan model -----
# model_string = "
# data{
#   int<lower=0> nRefBarcode ; // number of barcodes in ref allele
#   int<lower=0> nMutBarcode ;
#   int<lower=0> nDNAblocks ; // number of DNA replicates / blocks / samples / transfections
#   int<lower=0> nRNAblocks ;
#   int<lower=0> refDNAmat[nRefBarcode, nDNAblocks] ; // MPRA count matrix. Rows = barcodes, columns = samples/blocks
#   int<lower=0> refRNAmat[nRefBarcode, nRNAblocks] ;
#   int<lower=0> mutDNAmat[nMutBarcode, nDNAblocks] ;
#   int<lower=0> mutRNAmat[nMutBarcode, nRNAblocks] ;
#   real<lower=0> muRefRNAHyperParams[2, nRNAblocks] ; // gamma hyper-parameters on negative binomial parameters
#   real<lower=0> phiRefRNAHyperParams[2, nRNAblocks] ;
#   real<lower=0> muMutRNAHyperParams[2, nRNAblocks] ;
#   real<lower=0> phiMutRNAHyperParams[2, nRNAblocks] ;
#   real<lower=0> muDNAHyperParams[2, nDNAblocks] ;
#   real<lower=0> phiDNAHyperParams[2, nDNAblocks] ;
# }
# parameters{
#   real<lower=0> muRefDNA[nDNAblocks] ; //mean parameters for each block in each allele for each nucleic acid
#   real<lower=0> muRefRNA[nRNAblocks] ;
#   real<lower=0> muMutDNA[nDNAblocks] ;
#   real<lower=0> muMutRNA[nRNAblocks] ;
#   real<lower=0> phiRefDNA[nDNAblocks] ; // size parameters
#   real<lower=0> phiRefRNA[nRNAblocks] ;
#   real<lower=0> phiMutDNA[nDNAblocks] ;
#   real<lower=0> phiMutRNA[nRNAblocks] ;
# }
# model{
# 
# 
#   for (i in 1:nDNAblocks){
# 
#     // negative binomial parameters come from gamma hyper-priors
#     muRefDNA[i] ~ gamma(muDNAHyperParams[1, i], muDNAHyperParams[2, i]) ; 
#     phiRefDNA[i] ~ gamma(phiDNAHyperParams[1, i], phiDNAHyperParams[2, i]) ;
# 
#     // count data comes from the specified negative binomial
#     refDNAmat[,i] ~ neg_binomial_2(muRefDNA[i], phiRefDNA[i]) ;
#     mutDNAmat[,i] ~ neg_binomial_2(muMutDNA[i], phiMutDNA[i]) ;
#   }
# 
#   for (i in 1:nRNAblocks){
#     // negative binomial parameters come from gamma hyper-priors
#     muRefRNA[i] ~ gamma(muRefRNAHyperParams[1, i], muRefRNAHyperParams[2, i]) ;
#     muMutRNA[i] ~ gamma(muMutRNAHyperParams[1, i], muMutRNAHyperParams[2, i]) ;
#     
#     phiRefRNA[i] ~ gamma(phiRefRNAHyperParams[1, i], phiRefRNAHyperParams[2, i]) ;
#     phiMutRNA[i] ~ gamma(phiMutRNAHyperParams[1, i], phiMutRNAHyperParams[2, i]) ;
# 
#     // count data comes from the specified negative binomial
#     refRNAmat[,i] ~ neg_binomial_2(muRefRNA[i], phiRefRNA[i]) ;
#     mutRNAmat[,i] ~ neg_binomial_2(muMutRNA[i], phiMutRNA[i]) ;
#   }
# }
# "
# 
# model_object = stan_model(model_code = model_string)

#' Generate a distance matrix from a matrix of predictors
generateDistMat = function(predictors, log_distance = TRUE) {
  #predictors is a n x d data frame of predictors
  
  if (log_distance) {
    predictors %>% 
      as.data.frame() %>% 
      as.matrix %>% 
      dist %>% 
      as.matrix %>% 
      log1p
  } else {
    predictors %>% 
      as.data.frame() %>% 
      as.matrix %>% 
      dist %>% 
      as.matrix
  }
  
}

#' 
findWeights = function(i, distMat, minDistKernel, minNumContributing = 30, increaseFold = 1.333){
  # for the ith variant, produce a vector of weightings such that at minimum minNumContributing
  # variants meaningfully contribute (i.e. cumsum(sortedWeights)[30] <= .99)
  # arguments after the 2nd are heuristics that may need tuning
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix (precomputed for speed)
  distKernel = minDistKernel
  
  rawWeights = dnorm(distMat[i,-i], sd = distKernel)
  scaledWeights = rawWeights / sum(rawWeights)
  sorted = sort(scaledWeights, decreasing = TRUE)
  
  notEnoughContributing = cumsum(sorted[1:minNumContributing])[minNumContributing] > .99
  #allZero = all(rawWeights$x == 0)
  
  # if there aren't more than minNumContributing variants providing meaningful contribution to the prior
  if (is.na(notEnoughContributing) || notEnoughContributing) { 
    
    # iteratively increase the kernel bandwith until they do
    while (is.na(notEnoughContributing) || notEnoughContributing) {
      distKernel = distKernel * increaseFold
      rawWeights = dnorm(distMat[i,-i], sd = distKernel)
      scaledWeights = rawWeights / sum(rawWeights)
      sorted = sort(scaledWeights, decreasing = TRUE)
      
      notEnoughContributing = cumsum(sorted[1:minNumContributing])[minNumContributing] > .99
    }
  }
  
  scaledWeights
}

#' @title Bayesian analysis of MPRA data
#' @description Given MPRA data and a set of predictors, perform a Bayesian analysis of variants using an empirical prior
#' @param mpra_data a data frame of mpra data
#' @param predictors a matching data frame of annotations
#' @param marginal_prior logical indicating whether or not to disregard the functional predictors and use a marginal prior estimated from the entire assay
#' @details \code{mpra_data} must meet the following format conditions
#'   \enumerate { 
#'     \item one row per barcode 
#'     \item one column of variant IDs (e.g. rs IDs) 
#'     \item one column of alleles. These must be character strings of either "ref" or "mut" 
#'     \item one additional column for every sequencing sample 
#'     \item column names of plasmid library samples must contain "DNA" (e.g. "DNA_1", "DNA_2", ...) 
#'     \item column names of samples from transcription products must contain "RNA" (e.g. "RNA_1", "RNA_2", ...) }
bayesian_mpra_analyze = function(mpra_data, predictors, marginal_prior = FALSE) {
  
  distMat = generateDistMat(predictors)
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix
  minDistKernel = distMat[upper.tri(distMat)] %>% 
    unlist() %>% 
    sort() %>% #sort all observed distances
    .[. > 0] %>% 
    quantile(.001) # pick the .1th quantile. The only variants that will use this kernel will be in very densely populated regions of predictor space
  
}