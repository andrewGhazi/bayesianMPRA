# let's lay out the skeleton for the package mostly adapted over from
# stanModelInformativePrior and changed to be generalized functions where
# applicable

# library(tidyverse)
# library(magrittr)
# library(rstan)
# library(parallel)
# library(fitdistrplus)
# source('R/mledistModified.R')

# # ### Stan model -----
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
#     muMutDNA[i] ~ gamma(muDNAHyperParams[1, i], muDNAHyperParams[2, i]) ; 
#     phiMutDNA[i] ~ gamma(phiDNAHyperParams[1, i], phiDNAHyperParams[2, i]) ;
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
#'
#' Given an nxd matrix of variant annotations, produce an nxn distance matrix
#' describing the inter-variant distances in annotation space
#' 
#' @param predictors an n x d data frame of predictors
#' @param log_distance a logical indicating to use the log1p of the distances (TRUE) or the raw euclidean distances (FALSE)
#' 
#' @importFrom magrittr %>%
generate_dist_mat = function(predictors, log_distance = TRUE) {
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

#' Produce annotation-based weightings
#'
#' For the ith variant, adaptively produce a vector of weightings such that at minimum
#' min_num_contributing variants meaningfully contribute (i.e.
#' cumsum(sortedWeights)[30] <= .99) arguments after the 2nd are heuristics that
#' may need tuning
#'
#' @param i the index of the variant in the assay
#' @param dist_mat annotation-based distance matrix
#' @param min_dist_kernel the minimum kernel possible to use in the annotation
#'   space
#' @param min_num_contributing the minimum number of variants that must
#'   contribute to the ith variant's prior
#' @param increase_fold the multiplicative amount by which to increase the
#'   kernel in the case the current kernel doesn't allow at least
#'   \code{min_num_contributing} variants to contribute
#'   
#' @details The kernel is initialized at some small, pre-computed value then
#'   iteratively increased until there are "enough" variants contributing. This
#'   keeps the prior for one variant from being too strongly influenced by
#'   extremely close neighbors in annotation space.
#' @return a similarity-based vector of weights of the other n-1 variants
find_weights = function(i, dist_mat, min_dist_kernel, min_num_contributing = 30, increase_fold = 1.333){
  # for the ith variant, produce a vector of weightings such that at minimum min_num_contributing
  # variants meaningfully contribute (i.e. cumsum(sortedWeights)[30] <= .99)
  # arguments after the 2nd are heuristics that may need tuning
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix (precomputed for speed)
  dist_kernel = min_dist_kernel
  
  raw_weights = dnorm(dist_mat[i,-i], sd = dist_kernel)
  scaled_weights = raw_weights / sum(raw_weights)
  sorted = sort(scaled_weights, decreasing = TRUE)
  
  not_enough_contributing = cumsum(sorted[1:min_num_contributing])[min_num_contributing] > .99
  #allZero = all(raw_weights$x == 0)
  
  # if there aren't more than min_num_contributing variants providing meaningful contribution to the prior
  if (is.na(not_enough_contributing) || not_enough_contributing) { 
    
    # iteratively increase the kernel bandwith until they do
    while (is.na(not_enough_contributing) || not_enough_contributing) {
      dist_kernel = dist_kernel * increase_fold
      raw_weights = dnorm(dist_mat[i,-i], sd = dist_kernel)
      scaled_weights = raw_weights / sum(raw_weights)
      sorted = sort(scaled_weights, decreasing = TRUE)
      
      not_enough_contributing = cumsum(sorted[1:min_num_contributing])[min_num_contributing] > .99
    }
  }
  
  scaled_weights
}

#' Safely fit a negative binomial
#' 
#' Safely fit a negative binomial
#' 
#' @param count_vec a vector of counts
#' 
#' @return a 2-element list of a result or error (the other is NULL)
#' 
#' @importFrom purrr safely
#' @importFrom purrr set_names
#' @importFrom fitdistrplus fitdist
safely_fit_negbin = function(count_vec){
  #This doesn't quite suppress all error messages but it does output the right things
  #an error can still be printed when there's very low variability of low counts (e.g. c(1, rep(0, 14)))
  safely(fitdist, 
         otherwise = list(estimate = purrr::set_names(rep(NA, 2), nm = c('mu', 'size'))), 
         quiet = TRUE)(count_vec, 'nbinom')
}

#' Estimate Transfection Parameters
#' 
#' Estimate Transfection Parameters
#' 
#' @param count_dat 
estTransfectionParameters = function(count_dat){
  # count_dat - a tibble with a type (ref/mut) column and columns of observed MPRA counts in given transfections
  # uses fitdistrplus::fitdist because MASS::fitdistr was cracking wise at me
  # using a modified version of fitdistrplus::mledist because the regular version can't use non-integer weights
  
  count_dat %>% 
    gather(block, count, -type) %>% 
    group_by(type, block) %>% 
    summarise(MLEnegBin = list(safelyFitNegBin(count))) %>% 
    ungroup %>% 
    mutate(muEst = map_dbl(MLEnegBin, ~.x$result$estimate['mu']),
           sizeEst = map_dbl(MLEnegBin, ~.x$result$estimate['size'])) %>% 
    dplyr::select(-MLEnegBin)
}

#' @title Bayesian analysis of MPRA data
#' @description Given MPRA data and a set of predictors, perform a Bayesian analysis of variants using an empirical prior
#' @param mpra_data a data frame of mpra data
#' @param predictors a matching data frame of annotations
#' @param marginal_prior logical indicating whether or not to disregard the functional predictors and use a marginal prior estimated from the entire assay
#' @details \code{mpra_data} must meet the following format conditions:
#'   \enumerate { 
#'     \item one row per barcode 
#'     \item one column of variant IDs (e.g. rs IDs) 
#'     \item one column of alleles. These must be character strings of either "ref" or "mut" 
#'     \item one additional column for every sequencing sample 
#'     \item column names of plasmid library samples must contain "DNA" (e.g. "DNA_1", "DNA_2", ...) 
#'     \item column names of samples from transcription products must contain "RNA" (e.g. "RNA_1", "RNA_2", ...) }
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom purrr map
bayesian_mpra_analyze = function(mpra_data, predictors, marginal_prior = FALSE) {
  
  dist_mat = generate_dist_mat(predictors)
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix
  min_dist_kernel = dist_mat[upper.tri(dist_mat)] %>% 
    unlist() %>% 
    sort() %>% #sort all observed distances
    .[. > 0] %>% 
    quantile(.001) # pick the .1th quantile. The only variants that will use this kernel will be in very densely populated regions of predictor space
  
  mpra_data %<>% 
    group_by(snp_id) %>% 
    nest(.key = count_data) %>% 
    mutate(weights = map(1:nrow(.), ~find_weights(.x, dist_mat, min_dist_kernel)))
    
    
    
  mpra_data %<>% 
    mutate(weights = map(1:nrow(.), ~findWeights(.x, distMat, minDistKernel)))
  
}