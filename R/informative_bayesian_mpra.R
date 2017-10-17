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
# modelString = "
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
# model = stan_model(model_code = modelString)

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
  
  if (log_dist) {
    predictors %>% 
      dplyr::select(snp_id) %>% 
      as.data.frame() %>% 
      as.matrix %>% 
      scale %>% 
      dist %>% 
      as.matrix %>% 
      log1p
  } else {
    predictors %>% 
      dplyr::select(snp_id) %>% 
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
findWeights = function(i, dist_mat, min_dist_kernel, min_num_contributing = 30, increase_fold = 1.333){
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
fit_mu_gamma = function(weights, mu_estimates){
  fn_to_min = function(param_vec){-sum(weights * dgamma(mu_estimates, 
                                                        shape = param_vec[1], 
                                                        rate = param_vec[2], 
                                                        log = TRUE))}
  
  mean_est = mean(mu_estimates)
  var_est = var(mu_estimates)
  
  initial_guess = c(mean_est**2 / var_est, mean_est/var_est)
  
  optim_res = optim(initial_guess, 
                    fn_to_min, 
                    lower = c(1e-12), 
                    control = list(ndeps = c(1e-4, 1e-5)))
  
  # The ndeps control option is needed to scale the optimization proposals.
  # If the rate suggestions are too huge then fn_to_min throws out Inf which
  # breaks the optimizer
  
  if (optim_res$convergence != 0) {
    stop(paste0('problems with gamma fitting, convergence code: ', optim_res$convergence))
  }
  
  optim_res$par %>% 
    set_names(c('shape', 'rate'))
}

fit_size_gamma = function(weights, size_estimates){
  # different ndeps, no lower bound so the optimizer works
  fn_to_min = function(param_vec){-sum(weights * dgamma(size_estimates, 
                                                        shape = param_vec[1], 
                                                        rate = param_vec[2], 
                                                        log = TRUE))}
  mean_est = mean(size_estimates)
  var_est = var(size_estimates)
  
  initial_guess = c(mean_est**2 / var_est, mean_est/var_est)
  
  optim_res = optim(initial_guess, 
                    fn_to_min, 
                    control = list(ndeps = c(1e-4, 1e-4))) 
  
  if (optim_res$convergence != 0) {
    stop(paste0('problems with gamma fitting, convergence code: ', optim_res$convergence))
  }
  
  optim_res$par %>% 
    set_names(c('shape', 'rate'))
}

plus_or_homebrew_mu = function(weights, mu_estimates, initial_mu_guess){
  # Try to fit the gamma with fitdistrplus. If that doesn't work, try the
  # homebrew optimizer that calculates its own allele-block initial guess for the
  # individual allele-block
  
  res = try(fitdistMod(mu_estimates, 
                       'gamma', 
                       weights = weights, 
                       start = initial_mu_guess,
                       control = list(ndeps = c(1e-4, 1e-5)),
                       lower = 1e-12)$estimate,
            silent = TRUE)
  
  if (class(res) == 'try-error') {
    fit_mu_gamma(weights, mu_estimates)
  } else{
    res
  }
}

plus_or_homebrew_size = function(weights, size_estimates, initial_size_guess){
  #same as above just different ndeps
  
  res = try(fitdistMod(size_estimates, 
                       'gamma', 
                       weights = weights, 
                       start = initial_size_guess,
                       control = list(ndeps = c(1e-4, 1e-4)),
                       lower = 1e-12)$estimate,
            silent = TRUE)
  
  if (class(res) == 'try-error') {
    fit_size_gamma(weights, size_estimates)
  } else{
    res
  }
}

fit_gamma_priors = function(snp_id_num){ 
  snp_in_question = mpra_data[snp_id_num,]
  others = mpra_data[-snp_id_num,]
  
  
  data_for_estimates = others %>% 
    mutate(weight = snp_in_question$weights[[1]]) %>% 
    dplyr::select(snp_id, nb_params, weight) %>% 
    unnest %>% 
    na.omit %>% 
    filter(weight > 1e-4*sort(snp_in_question$weights[[1]], decreasing = TRUE)[30], # This is necessary for speed)
           grepl('RNA', block)) # Only fit conditional prior on RNA. DNA counts use a marginal prior
  
  mu_dat = data_for_estimates$nb_mu_est
  initial_mu_guess = list(shape = mean(mu_dat)**2 / var(mu_dat), 
                          rate = mean(mu_dat) / var(mu_dat))
  
  size_dat = data_for_estimates$nb_size_est[data_for_estimates$nb_size_est < quantile(data_for_estimates$nb_size_est, .99)]
  initial_size_guess = list(shape = mean(size_dat)**2 / var(size_dat), 
                            rate = mean(size_dat) / var(size_dat))
  
  data_for_estimates %>% 
    group_by(allele, block) %>% 
    summarise(mu_gamma_priors = list(plus_or_homebrew_mu(weight, 
                                                         nb_mu_est, 
                                                         initial_mu_guess)),
              size_gamma_priors = list(plus_or_homebrew_size(weight[nb_size_est < 1e4], # Have to impose a size restriction / magic number
                                                             nb_size_est[nb_size_est < 1e4], 
                                                             initial_size_guess))) %>% 
    ungroup
  
  # Estimates of the size parameter can be unstable (because variance = mu + mu^2
  # / size so size --> Inf as var --> mu) so we cut out those that are above the
  # 99th quantile. These are variants that are essentially poisson in their
  # counts so we're slightly biasing our result to show HIGHER variance.
}

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
#' @importFrom magrittr %<>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom purrr reduce
#' @importFrom rstan sampling
#' @importFrom coda mcmc
#' @importFrom coda HPDinterval
run_sampler = function(snp_data, marg_dna_priors, save_nonfunctional, out_dir, num_chain = 3, num_iter = 3334, num_warmup = 500){
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
  
  mu_DNA_hyper_params = margDNAPrior %>% 
    filter(grepl('mu', negBinParam)) %>% 
    dplyr::select(alpha_est, beta_est) %>% # OMG three pipes lining up
    as.matrix %>% 
    t
  
  phi_DNA_hyper_params = margDNAPrior %>% 
    filter(grepl('size', negBinParam)) %>% 
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
  
  sampler_res = sampling(object = model,
                         data = data_list,
                         chains = num_chain,
                         iter = num_iter,
                         warmup = num_warmup,
                         thin = 1,
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
         file = paste0(out_dir, snp_data$snp_id, '.RData'))
  }
}


#' @title Bayesian analysis of MPRA data
#' @description Given MPRA data and a set of predictors, perform a Bayesian analysis of variants using an empirical prior
#' @param mpra_data a data frame of mpra data
#' @param predictors a matching data frame of annotations
#' @param out_dir a directory that you want the outputs written to. Make sure it ends with a forward slash.
#' @param save_nonfunctional logical indicating whether to save the sampler results of non-functional variants. 
#' @param marginal_prior logical indicating whether or not to disregard the functional predictors and use a marginal prior estimated from the entire assay
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
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom purrr map
#' @importFrom purrr nest
#' @export
bayesian_mpra_analyze = function(mpra_data, predictors, use_marg_prior = FALSE, out_dir, save_nonfunctional = FALSE, normalization_method = 'quantile_normalization', num_cores = 1) {
  # mpra_data is a data frame with columns like so:
  # one column called snp_id
  
  if (!use_marg_prior) {
    # Initialize the kernel at some small value based on the typical distances in the input distance matrix
    ordered_preds = mpra_data %>% 
      dplyr::select(snp_id) %>% 
      unique %>% 
      left_join(predictors, by = 'snp_id') %>% 
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
      mutate(weights = mclapply(1:nrow(.), findWeights, dist_mat = dist_mat, min_dist_kernel = min_dist_kernel, mc.cores = num_cores))
    
    print('Computing neg_bin parameters by sample and allele')
    mpra_data %<>%
      mutate(nb_params = mclapply(count_data, est_sample_params, mc.cores = num_cores))
    
    print('Fitting annotation-based gamma prior...')
    mpra_data %<>% 
      mutate(rna_gamma_params = mclapply(1:n(), fit_gamma_priors, mc.cores = num_cores))
  } else {
    # assign marginal priors to the RNA samples
    stop('Marginal DNA priors are not yet implemented.')
  }
  
  print('Fitting marginal prior to DNA samples...')
  marg_dna_prior = fit_DNA_prior(mpra_data)
  
  print('Running MCMC samplers...')
  mpra_data %>% 
    group_by(snp_id) %>%
    nest %>%
    mutate(sampler_result = mclapply(data, run_sampler, 
                                     marg_dna_prior = marg_dna_prior, 
                                     save_nonfunctional = save_nonfunctional,
                                     out_dir = out_dir,
                                     norm_method = normalization_method,
                                     mc.cores = num_cores))
  
}