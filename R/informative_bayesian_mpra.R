# let's lay out the skeleton for the package mostly adapted over from
# stanModelInformativePrior and changed to be generalized functions where
# applicable

# library(tidyverse)
# library(magrittr)
# library(rstan)
# library(parallel)
# library(fitdistrplus)
# source('analysis/mledistModified.R')
# 
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
#' @details the predictors need to be ordered by snp_id in the same order as they appear in mpra_data
#' 
#' @importFrom magrittr %>%
generate_dist_mat = function(predictors, log_distance = TRUE) {
  #predictors is a n x d data frame of predictors
  
  if (log_distance) {
    predictors %>% 
      as.data.frame() %>% 
      as.matrix %>% 
      scale %>% 
      dist %>% 
      as.matrix %>% 
      log1p
  } else {
    predictors %>% 
      as.data.frame() %>% 
      as.matrix %>% 
      scale %>% 
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
#' 
find_weights = function(i, dist_mat, min_dist_kernel, min_num_contributing = 30, increase_fold = 1.333){
  # for the ith variant, produce a vector of weightings such that at minimum min_num_contributing
  # variants meaningfully contribute (i.e. cumsum(sortedWeights)[30] <= .99)
  # arguments after the 2nd are heuristics that may need tuning
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix (precomputed for speed)
  dist_kernel = min_dist_kernel
  
  raw_weights = dnorm(dist_mat[i,-i], sd = dist_kernel)
  scaled_weights = raw_weights / sum(raw_weights, na.rm = TRUE)
  sorted = sort(scaled_weights, decreasing = TRUE)
  
  not_enough_contributing = cumsum(sorted[1:min_num_contributing])[min_num_contributing] > .99
  #allZero = all(raw_weights$x == 0)
  
  # if there aren't more than min_num_contributing variants providing meaningful contribution to the prior
  if (is.na(not_enough_contributing) || not_enough_contributing) { 
    
    # iteratively increase the kernel bandwith until they do
    while (is.na(not_enough_contributing) || not_enough_contributing) {
      dist_kernel = dist_kernel * increase_fold
      raw_weights = dnorm(dist_mat[i,-i], sd = dist_kernel)
      scaled_weights = raw_weights / sum(raw_weights, na.rm = TRUE)
      sorted = sort(scaled_weights, decreasing = TRUE)
      
      not_enough_contributing = cumsum(sorted[1:min_num_contributing])[min_num_contributing] > .99
    }
  }
  
  scaled_weights
}

### Fit NegBin estimates to each variant -----
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
safelyFitNegBin = function(count_vec){
  #This doesn't quite suppress all error messages but it does output the right things
  #an error can still be printed when there's very low variability of low counts (e.g. c(1, rep(0, 14)))
  safely(fitdist, 
         otherwise = list(estimate = purrr::set_names(rep(NA, 2), nm = c('mu', 'size'))), 
         quiet = TRUE)(count_vec, 'nbinom')
}


#' Estimate Transfection Parameters
#' 
#' Estimate Transfection Parameters for a single allele
#' 
#' @param count_dat 
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#' @importFrom tidyr gather
#' 
est_sample_params = function(count_dat){
  # count_dat - a tibble with a allele (ref/mut) column and columns of observed MPRA counts in given transfections
  # uses fitdistrplus::fitdist because MASS::fitdistr was cracking wise at me
  # using a modified version of fitdistrplus::mledist because the regular version can't use non-integer weights (and it was cracking wise)
  
  count_dat %>% 
    gather(block, count, -allele) %>% 
    group_by(allele, block) %>% 
    summarise(mle_neg_bin = list(safelyFitNegBin(count))) %>% 
    ungroup %>% 
    mutate(nb_mu_est = map_dbl(mle_neg_bin, ~.x$result$estimate['mu']),
           nb_size_est = map_dbl(mle_neg_bin, ~.x$result$estimate['size'])) %>% 
    dplyr::select(-mle_neg_bin)
}

### Fit weighted gamma hyperprior on negBin parameters to each variant ------- 

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
fit_gamma_priors = function(snp_id_num, mpra_data, marg_rna_gamma_prior){ 
  snp_in_question = mpra_data[snp_id_num,]
  others = mpra_data[-snp_id_num,]
  
  data_for_estimates = others %>% 
    mutate(weight = snp_in_question$weights[[1]]) %>% 
    dplyr::select(snp_id, nb_params, weight) %>% 
    unnest %>% 
    na.omit
  
  gamma_priors = data_for_estimates %>% 
    group_by(allele, block) %>% 
    summarise(mu_gamma_priors = list(try(fitdist(nb_mu_est, 'gamma', lower = c(0,0), weights = weight))),
              size_gamma_priors = list(try(fitdist(nb_size_est, 'gamma', lower = c(0,0), weights = weight)))) %>% 
    left_join(marg_rna_gamma_prior)
  
  ## If the gamma fitting failed, find where
  failed_mu = map_lgl(gamma_priors$mu_gamma_priors, ~class(.x) == 'try-error')
  failed_size = map_lgl(gamma_priors$size_gamma_priors, ~class(.x) == 'try-error')
  
  # Replace failed prior estimates with marginal estimates
  gamma_priors$mu_gamma_priors[failed_mu] = gamma_priors$mu_marg_gamma_prior[failed_mu]
  gamma_priors$size_gamma_priors[failed_size] = gamma_priors$size_marg_gamma_prior[failed_size]
  
  gamma_priors %>% 
    mutate(mu_gamma_priors = map(mu_gamma_priors, ~.x$estimate),
           size_gamma_priors = map(size_gamma_priors, ~.x$estimate)) %>% 
    select(allele, block, mu_gamma_priors, size_gamma_priors)
  
}

#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @importFrom dplyr filter
#' @importFrom tidyr gather
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom purrr map_dbl
#' @importFrom fitdistrplus fitdist
fit_DNA_prior = function(mpra_data){
  mpra_data %>% 
    dplyr::select(snp_id, nb_params) %>% 
    unnest %>% 
    filter(grepl('DNA', block)) %>%
    gather(nb_param, nb_param_val, -(snp_id:block)) %>% 
    group_by(block, nb_param) %>% 
    summarise(marg_gamma_estimate = list(fitdist(nb_param_val, 'gamma'))) %>% 
    ungroup %>% 
    mutate(alpha_est = map_dbl(marg_gamma_estimate, ~.x$estimate[1]),
           beta_est = map_dbl(marg_gamma_estimate, ~.x$estimate[2]))
}


#' quantile normalize parameters
#' 
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom dplyr select
#' @importFrom dplyr contains
#' @importFrom magrittr %>%
#' @importFrom tibble as.tibble
#' @importFrom purrr set_names
quant_norm_parameter = function(param_name, sample_df){
  sample_df %>%
    dplyr::select(contains(param_name)) %>%
    as.matrix %>%
    normalize.quantiles %>%
    rowMeans %>% # Still not sure that taking the mean across samples is the right thing to do
    as.tibble %>%
    set_names(param_name)
}

#' Compute transcriptional shift samples
#' 
#' Compute transcriptional shift samples from sampler output
#' 
#' @importFrom magrittr %>%
#' @importFrom rstan extract
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom purrr set_names
#' @importFrom dplyr bind_cols
#' @importFrom dplyr transmute
get_snp_TS_samples = function(sampler_res){
  sample_df = sampler_res %>%
    rstan::extract() %>%
    map(as.tibble) %>%
    map2(names(.),
         .,
         ~set_names(.y, paste0(.x, '_', names(.y)))) %>%
    bind_cols
  
  c('muMutRNA',
    'muMutDNA',
    'muRefRNA',
    'muRefDNA') %>%
    map(quant_norm_parameter, sample_df = sample_df) %>%
    bind_cols %>%
    transmute(transcriptional_shift = log(muMutRNA / muMutDNA) - log(muRefRNA / muRefDNA))
}

#' Run the Stan Sampler
#' 
#' Run the Stan sampler for one variant
#' 
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr contains
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom purrr reduce
#' @importFrom rstan sampling
#' @importFrom coda mcmc
#' @importFrom coda HPDinterval
#' @importFrom preprocessCore normalize.quantiles
run_sampler = function(snp_data, marg_dna_prior, save_nonfunctional, out_dir, num_chain = 3, num_iter = 3334, num_warmup = 500, norm_method = 'quantile_normalization', object){
  # snp_data - a data_frame with one row containing a column called count_data and another called rna_gamma_params and a column called snp_id
  
  # Given a matrix of counts (rows = barcodes, columns = samples) and a
  # data_frame of by-allele-RNA-sample gamma hyperpriors, run the above Stan
  # model
  
  count_data = snp_data$count_data[[1]]
  RNA_gamma_params = snp_data$rna_gamma_params[[1]]
  
  # Prepare count data matrices
  ref_DNA_mat = count_data %>% 
    filter(allele == 'ref') %>% 
    dplyr::select(-allele) %>% 
    dplyr::select(contains('DNA')) %>% 
    as.matrix
  
  ref_RNA_mat = count_data %>% 
    filter(allele == 'ref') %>% 
    dplyr::select(-allele) %>% 
    dplyr::select(contains('RNA')) %>% 
    as.matrix
  
  mut_DNA_mat  = count_data %>% 
    filter(allele == 'mut') %>% 
    dplyr::select(-allele) %>% 
    dplyr::select(contains('DNA')) %>% 
    as.matrix
  
  mut_RNA_mat = count_data %>% 
    filter(allele == 'mut') %>% 
    dplyr::select(-allele) %>% 
    dplyr::select(contains('RNA')) %>% 
    as.matrix
  
  # Prepare Gamma hyper-prior matrices
  mu_ref_RNA_hyper_params = RNA_gamma_params %>% 
    filter(allele == 'ref') %>% 
    pull(mu_gamma_priors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  mu_mut_RNA_hyper_params = RNA_gamma_params %>% 
    filter(allele == 'mut') %>% 
    pull(mu_gamma_priors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  phi_ref_RNA_hyper_params = RNA_gamma_params %>% 
    filter(allele == 'ref') %>% 
    pull(size_gamma_priors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  phi_mut_RNA_hyper_params = RNA_gamma_params %>% 
    filter(allele == 'mut') %>% 
    pull(size_gamma_priors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  mu_DNA_hyper_params = marg_dna_prior %>% 
    filter(grepl('mu', nb_param)) %>% 
    dplyr::select(alpha_est, beta_est) %>% # OMG three pipes lining up
    as.matrix %>% 
    t
  
  phi_DNA_hyper_params = marg_dna_prior %>% 
    filter(grepl('size', nb_param)) %>% 
    dplyr::select(alpha_est, beta_est) %>% 
    as.matrix %>% 
    t
  
  # create input data list
  data_list = list(nRefBarcode = nrow(ref_DNA_mat),
                   nMutBarcode = nrow(mut_DNA_mat), 
                   nDNAblocks = ncol(ref_DNA_mat), 
                   nRNAblocks = ncol(ref_RNA_mat),
                   refDNAmat = ref_DNA_mat,
                   refRNAmat = ref_RNA_mat,
                   mutDNAmat = mut_DNA_mat,
                   mutRNAmat = mut_RNA_mat,
                   muRefRNAHyperParams = mu_ref_RNA_hyper_params,
                   phiRefRNAHyperParams = phi_ref_RNA_hyper_params,
                   muMutRNAHyperParams = mu_mut_RNA_hyper_params,
                   phiMutRNAHyperParams = phi_mut_RNA_hyper_params,
                   muDNAHyperParams = mu_DNA_hyper_params,
                   phiDNAHyperParams = phi_DNA_hyper_params)
  
  sampler_res = sampling(object = object,
                         data = data_list,
                         chains = num_chain,
                         iter = num_iter,
                         warmup = num_warmup,
                         thin = 1,
                         cores = 1,
                         verbose = FALSE) #friggin stan still verbose af
  
  if (norm_method == 'quantile_normalization') {
    ts_samples = get_snp_TS_samples(sampler_res)
    
    ts_hdi = ts_samples %>% mcmc %>% HPDinterval(prob = .95)
    
    functional_variant = !between(0, ts_hdi[1], ts_hdi[2])
    
  } else {
    stop('Normalization methods other than quantile normalization are not yet implemented yet.')
  }
  
  if (functional_variant) {
    mean_transcriptional_shift = mean(ts_samples)
    
    save(sampler_res,
         ts_samples,
         ts_hdi,
         mean_transcriptional_shift,
         file = paste0(out_dir, gsub(' ', '_', gsub('/', '-', snp_data$snp_id)), '.RData'))
    return('functional')
  } else {
    return('non-functional')
  }
}


#' @title Bayesian analysis of MPRA data
#' @description Given MPRA data and a set of predictors, perform a Bayesian analysis of variants using an empirical prior
#' @param mpra_data a data frame of mpra data
#' @param predictors a matching data frame of annotations
#' @param use_marg_prior logical indicating whether or not to disregard the functional predictors and use a marginal prior estimated from the entire assay
#' @param out_dir a directory that you want the outputs written to. Make sure it ends with a forward slash.
#' @param mpra_model_object a stan model object compiled with rstan::stan_model(model_code = mpra_model_string). mpra_model_string is a data object bundled with the package.
#' @param save_nonfunctional logical indicating whether to save the sampler results of non-functional variants. 
#' @param normalization_method character vector indicating which method to use for aggregating information across samples. Must be either 'quantile_normalization' or 'depth_normalization'
#' @param num_cores integer indicating how many cores to use for parallelization. Currently the analysis takes ~15s per variant on a first-gen i7 CPU, so setting this as high as possible is recommended as long as you have plenty of RAM.
#' @details \code{mpra_data} must meet the following format conditions: \enumerate{ 
#'     \item one row per barcode 
#'     \item one column of variant IDs (e.g. rs IDs) 
#'     \item one column of alleles called `allele`. These must be character strings of either "ref" or "mut" 
#'     \item one additional column for every transfection/physical sample
#'     \item column names of plasmid library samples must contain "DNA" (e.g. "DNA_1", "DNA_2", ...) 
#'     \item column names of samples from transcription products must contain "RNA" (e.g. "RNA_1", "RNA_2", ...) 
#'     }
#'     
#'     \code{save_nonfunctional} defaults to \code{FALSE} as doing so can consume a large amount of storage space
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr matches
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom parallel mclapply
#' @importFrom purrr map
#' @importFrom purrr map2
#' @importFrom tidyr nest
#' @importFrom rstan stan_model
#' @export
bayesian_mpra_analyze = function(mpra_data, 
                                 predictors, 
                                 use_marg_prior = FALSE, 
                                 out_dir, 
                                 mpra_model_object,
                                 save_nonfunctional = FALSE,
                                 normalization_method = 'quantile_normalization', 
                                 num_cores = 1) {
  # mpra_data is a data frame with columns like so:
  # one column called snp_id
  
  if (missing(mpra_data)) {
    stop('mpra_data is missing: You must provide MPRA counts to do Bayesian MPRA analysis ಠ_ಠ')
  }
  
  if (missing(predictors)) {
    use_marg_prior = TRUE
    warning('predictors argument is empty. Argument use_marg_prior has been set to TRUE')
  }
  
  if (normalization_method != 'quantile_normalization') {
    stop('Normalization methods other than quantile normalization are not yet implemented.')
  }
  
  if (missing(out_dir)) {
    stop('You must specify an output directory out_dir')
  }
  
  if (use_marg_prior) {
    #stop('The use of a marginal prior is not yet implemented.')
    
    print('Organizing count data...')
    mpra_data %<>% 
      dplyr::select(matches('snp_id|RNA|DNA|allele')) %>% 
      group_by(snp_id) %>% 
      nest(.key = count_data) 
    
    print('Computing neg_bin parameters by sample and allele...')
    mpra_data %<>%
      mutate(nb_params = mclapply(count_data, est_sample_params, mc.cores = num_cores)) # ~11 / second / core
    
    print('Fitting marginal gamma priors on RNA samples...')
    marg_rna_gamma_priors = mpra_data %>%
      select(snp_id, nb_params) %>% 
      unnest() %>%
      na.omit %>% 
      filter(grepl('RNA', block)) %>% # The marginal DNA prior shouldn't care about ref/alt and it gets attached later
      group_by(block, allele) %>% 
      summarise(gamma_mu_prior = list(try(fitdist(data = nb_mu_est, distr = 'gamma', lower = 0))),
                gamma_size_prior = list(try(fitdist(data = nb_size_est, distr = 'gamma', lower = 0))),
                mu_gamma_priors = map(gamma_mu_prior, ~.x$estimate), # This and the next line just reformat things to work with the existing code
                size_gamma_priors = map(gamma_size_prior, ~.x$estimate)) %>% 
      select(-gamma_mu_prior, -gamma_size_prior) %>% 
      ungroup
    
    mpra_data %<>% 
      mutate(rna_gamma_params = list(marg_rna_gamma_priors))
    
  } else {
    # Initialize the kernel at some small value based on the typical distances in the input distance matrix
    ordered_preds = mpra_data %>% 
      dplyr::select(snp_id) %>% 
      unique %>% 
      left_join(predictors, by = 'snp_id') %>% 
      na.omit # drop snp_id's without annotations
    
    mpra_data %<>% 
      filter(snp_id %in% ordered_preds$snp_id) # drop snp_id's without annotations
    
    ordered_preds %<>%
      dplyr::select(-snp_id)
    
    print('Computing distance matrix...')
    dist_mat = generate_dist_mat(ordered_preds)
    
    min_dist_kernel = dist_mat[upper.tri(dist_mat)] %>% 
      unlist() %>% 
      sort() %>% #sort all observed distances
      .[. > 0] %>% 
      quantile(.001) # pick the .1th quantile. The only variants that will use this kernel will be in very densely populated regions of predictor space
    
    print('Organizing count data...')
    mpra_data %<>% 
      dplyr::select(matches('snp_id|RNA|DNA|allele')) %>% 
      group_by(snp_id) %>% 
      nest(.key = count_data) 
    
    print('Evaluating predictor based weights...')
    mpra_data %<>%
      mutate(weights = mclapply(1:nrow(.), find_weights, dist_mat = dist_mat, min_dist_kernel = min_dist_kernel, mc.cores = num_cores))
    
    print('Computing neg_bin parameters by sample and allele...')
    mpra_data %<>%
      mutate(nb_params = mclapply(count_data, est_sample_params, mc.cores = num_cores))
    
    marg_rna_gamma_prior = mpra_data %>% 
      dplyr::select(nb_params) %>% 
      unnest %>% 
      na.omit %>% 
      filter(grepl('RNA', block)) %>% 
      group_by(allele, block) %>% 
      summarise(mu_marg_gamma_prior = list(fitdist(data = nb_mu_est,
                                              distr = 'gamma',
                                              lower = c(0,0))),
                size_marg_gamma_prior = list(fitdist(data = nb_size_est,
                                                distr = 'gamma',
                                                lower = c(0,0))),
                mu_prior_shape = map_dbl(mu_marg_gamma_prior, ~.x$estimate[1]),
                mu_prior_rate = map_dbl(mu_marg_gamma_prior, ~.x$estimate[2]),
                size_prior_shape = map_dbl(size_marg_gamma_prior, ~.x$estimate[1]),
                size_prior_rate = map_dbl(size_marg_gamma_prior, ~.x$estimate[2]))
    
    print('Fitting annotation-based gamma prior...')
    mpra_data %<>%
      mutate(rna_gamma_params = mclapply(1:n(), 
                                         fit_gamma_priors, 
                                         mpra_data = mpra_data,
                                         marg_rna_gamma_prior = marg_rna_gamma_prior,
                                         mc.cores = num_cores))
  }
  
  print('Fitting marginal prior to DNA samples...')
  marg_dna_prior = fit_DNA_prior(mpra_data)
  
  print('Running MCMC samplers...')
  
  #model_object2 = stan_model(model_code = model_string)
  mpra_data %>% 
    group_by(snp_id) %>%
    nest %>%
    mutate(data = map2(snp_id, data, ~mutate(.y, snp_id = .x)),
           sampler_result = unlist(mclapply(data, 
                                            run_sampler, 
                                            marg_dna_prior = marg_dna_prior, 
                                            save_nonfunctional = save_nonfunctional,
                                            out_dir = out_dir,
                                            norm_method = normalization_method,
                                            object = mpra_model_object,
                                            mc.cores = num_cores)))
  
}
