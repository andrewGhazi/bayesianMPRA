library(tidyverse)
library(magrittr)
library(stringr)
library(rstan)
library(coda)

n_cores = 5

dir = "/mnt/labhome/andrew/MPRA/paper_data/"

UMPRA = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                   delim = "\t",
                   col_names = TRUE,
                   col_types = cols(chr = "c")) %>% 
  select(snp_id = construct, allele = type, matches('DNA|CTRL')) %>% 
  mutate(allele = gsub('Mut', 'alt', gsub('Ref', 'ref', allele)),
         bc_id = 1:n())

## Calculate depth factors and well-represented barcodes

depth_factors = UMPRA %>% 
  select(matches('[DR]NA')) %>% 
  map_df(~sum(.x) / 1e6) %>% 
  gather(sample, depth_factor)

depth_norm_dna = UMPRA %>% 
  gather(sample, count, matches('[DR]NA')) %>% 
  filter(grepl('DNA', sample)) %>% 
  left_join(depth_factors, by = 'sample') %>% 
  mutate(depth_norm_count = count / depth_factor) %>% 
  group_by(bc_id) %>% 
  summarise(mean_dn_dna = mean(depth_norm_count, na.rm = TRUE))

well_represented = depth_norm_dna %>% 
  filter(mean_dn_dna >= .13)

depth_norm_dna %>% 
  ggplot(aes(mean_dn_dna)) + 
  geom_histogram(aes(y = ..density..)) + 
  geom_density() + 
  geom_vline(xintercept = .13) + 
  geom_rug(alpha = .05) + 
  scale_x_log10()
  

## Estimate priors

fit_nb = function(counts){
  counts = counts$count
  fn_to_min = function(param_vec){
    # param_vec[1] nb mean
    # param_vec[2] nb size
    -sum(dnbinom(counts,
                 mu = param_vec[1],
                 size = param_vec[2],
                 log = TRUE))
  }
  
  stats::nlminb(start = c(100, 1),
                objective = fn_to_min,
                lower = rep(.Machine$double.xmin, 2))
}

fit_gamma = function(param_estimates){
  
  fn_to_min = function(ab_vec){
    -sum(dgamma(param_estimates,
                shape = ab_vec[1],
                rate = ab_vec[2],
                log = TRUE))
  }
  
  stats::nlminb(start = c(1,1),
                objective = fn_to_min,
                lower = rep(.Machine$double.xmin, 2))
}


nb_param_estimates = UMPRA %>% 
  filter(bc_id %in% well_represented$bc_id) %>% 
  gather(sample, count, matches('[DR]NA')) %>% 
  group_by(snp_id, allele, sample) %>% 
  nest %>%
  mutate(count_mean = map_dbl(data, ~mean(.x$count)),
         nb_fit = mclapply(data, fit_nb, mc.cores = n_cores),
         converged = map_lgl(nb_fit, ~.x$convergence == 0)) %>% 
  filter(converged) %>%
  left_join(depth_factors, by = 'sample') %>%
  mutate(depth_adj_mean = count_mean / depth_factor,
         depth_adj_mu_est = map2_dbl(nb_fit, depth_factor, ~.x$par[1] / .y),
         phi_est = map_dbl(nb_fit, ~.x$par[2])) 

marg_priors = nb_param_estimates %>%
  mutate(acid_type = factor(str_extract(sample, 'DNA|RNA'))) %>%
  group_by(allele, acid_type) %>%
  summarise(phi_gamma_prior = list(fit_gamma(phi_est)),
            mu_gamma_prior = list(fit_gamma(depth_adj_mu_est))) %>%
  ungroup %>%
  gather(prior_type, gamma_fit, matches('gamma')) %>%
  mutate(alpha_est = map_dbl(gamma_fit, ~.x$par[1]),
         beta_est = map_dbl(gamma_fit, ~.x$par[2])) %>%
  arrange(desc(allele)) 

#### Run the model on each allele
# only return the HPD on TS for each allele


my_HPD <- function(obj, prob = 0.95, ...) {
  dimnames(obj) = NULL # Stan outputs the iterations as only one dimname which makes as.matrix() fail
  obj <- as.matrix(obj)
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals)) stop("obj must have nsamp > 1")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                2, which.min)
  ans <- cbind(vals[cbind(inds, 1:npar)],
               vals[cbind(inds + gap, 1:npar)])
  dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
  attr(ans, "Probability") <- gap/nsamp
  ans
}



run_full_model = function(snp_id, 
                          snp_dat,
                          depth_factors,
                          barcode_prior,
                          marg_dna_prior,
                          marg_rna_prior, 
                          tot_samp,
                          n_cores,
                          n_warmup){
  
  
  # Prepare the stan model input data list
  data_list = list(n_rna_samples = snp_dat %>% names %>% grepl('RNA', .) %>% sum,
                   n_dna_samples = snp_dat %>% names %>% grepl('DNA', .) %>% sum,
                   n_ref = snp_dat %>% filter(tolower(allele) == 'ref') %>% nrow,
                   n_alt = snp_dat %>% filter(tolower(allele) != 'ref') %>% nrow,
                   ref_dna_counts = snp_dat %>% filter(tolower(allele) == 'ref') %>% select(matches('DNA')) %>% as.matrix,
                   alt_dna_counts = snp_dat %>% filter(tolower(allele) != 'ref') %>% select(matches('DNA')) %>% as.matrix,
                   ref_rna_counts = snp_dat %>% filter(tolower(allele) == 'ref') %>% select(matches('RNA')) %>% as.matrix,
                   alt_rna_counts = snp_dat %>% filter(tolower(allele) != 'ref') %>% select(matches('RNA')) %>% as.matrix,
                   rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
                   dna_depths = depth_factors %>% filter(grepl('DNA', sample)) %>% pull(depth_factor),
                   dna_m_a = dna_gamma_estimates$par[1],
                   dna_m_b = dna_gamma_estimates$par[2],
                   dna_p_a = marg_dna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(alpha_est), 
                   dna_p_b = marg_dna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(beta_est),
                   rna_m_a = marg_rna_prior %>% filter(prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
                   rna_m_b = marg_rna_prior %>% filter(prior_type == 'mu_gamma_prior') %>% pull(beta_est),
                   rna_p_a = marg_rna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(alpha_est), 
                   rna_p_b = marg_rna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(beta_est))
  
  # Calculate sampler parameters
  n_samp_per_core = tot_samp / n_cores 
  n_iter = n_samp_per_core + n_warmup
  
  # Run the sampler
  
  vb_res = vb()
  # sampler_res = sampling(bc_object,
  #                        data = data_list,
  #                        chains = n_cores, 
  #                        iter = n_iter, 
  #                        warmup = n_warmup,
  #                        cores = n_cores)
  # 
  # Save the sampler results
  save(vb_res,
       file = paste0('/mnt/bigData2/andrew/analysis_outputs/ulirsch_fully_parameterised/full_model_results/', snp_id, '.RData'))
  
  # Return a 95% HDI and mean and 99% HDI
  ts_samples = sampler_res %>% 
    rstan::extract(pars = 'transcription_shift') %>% 
    .[['transcription_shift']] %>% 
    coda::mcmc()
  
  hdi_95 = ts_samples %>% 
    my_HPD(prob = .95)
  
  hdi_99 = ts_samples %>% 
    my_HPD(prob = .99)
}