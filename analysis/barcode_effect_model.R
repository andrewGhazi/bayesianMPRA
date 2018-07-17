# I owe you the barcode normalization
#
# It goes like this .  Take n values.  Take the geometric mean of these values .
# Take the ith value .  Divide the ith value by the geometric mean .
# Reciprocate that .
#
# If you use the reciprocate value above as a normalizing scale factor you can
# remove the barcode effects .
#
# For the n values I have in mind the library normalized median values of each
# barcode in the collection
#
# So for a given barcode , 1)  take the depth normalized values in each
# replicate . 2) take the median of these values 3) apply the geometric mean
# based procedure above 4). Take the scale factors for barcodes that result and
# stick them into your modeling scheme
#
#
# Please let me know if you have questions or what you think
#
#
# Chad

library(tidyverse)
library(magrittr)
library(rstan)

# load("/mnt/bigData2/andrew/analysis_data/testing_dat.RData")
# 
# sample_depths = mpra_data %>%
#   unnest %>%
#   select(-snp_id, -allele) %>%
#   gather(sample, count) %>%
#   group_by(sample) %>%
#   summarise(depth = sum(count))
# 
# snp_dat = mpra_data$count_data[[1]] %>%
#   mutate(bc_id = 1:n())
# 
# # For one barcode take the depth normalized values in each replicate .
# 
# dnv = snp_dat %>%
#   gather(sample, count, -allele, -bc_id) %>%
#   left_join(sample_depths, by = 'sample') %>%
#   mutate(depth_norm_count = 1e6 * count / depth)
# 
# well_represented = dnv %>%
#   filter(grepl('DNA', sample)) %>%
#   group_by(allele, bc_id) %>%
#   summarise(mean_depth_norm = mean(depth_norm_count)) %>%
#   ungroup %>%
#   filter(mean_depth_norm > 10)
# 
# dnv %<>% filter(bc_id %in% well_represented$bc_id)
# 
# # 2) take the median of these values
# medians = dnv %>%
#   mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>%
#   group_by(bc_id, samp_type) %>%
#   summarise(med_dnv = median(depth_norm_count)) %>%
#   ungroup %>%
#   left_join(unique(select(dnv, bc_id, allele)),
#             by = 'bc_id')
# 
# # 3) apply the geometric mean based procedure above 
# #bc_medians = 
# 
# # https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
# gm_mean = function(x, na.rm=TRUE){ 
#   exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
# }
# 
# geo_means = medians %>%
#   group_by(allele, samp_type) %>%
#   summarise(geo_mean = gm_mean(med_dnv))
# 
# bc_norm_factors = medians %>%
#   left_join(geo_means, by = c('allele', 'samp_type')) %>%
#   mutate(bc_norm_factor = (med_dnv / geo_mean)^(-1)) %>%
#   select(bc_id, samp_type, bc_norm_factor) %>%
#   mutate(samp_type = paste0(samp_type, '_norm_factor')) %>%
#   spread(samp_type, bc_norm_factor)
# 
# inputs = snp_dat %>%
#   left_join(bc_norm_factors,
#             by = 'bc_id') %>%
#   select(allele, bc_id, DNA_norm_factor, RNA_norm_factor, everything())
# 
# ## prior estimation copied from neg_bin_regression.R ----
# 
# load('/mnt/bigData2/andrew/analysis_data/testing_dat_full.RData')
# mpra_data = ulirschCounts %>%
#   group_by(snp_id) %>%
#   nest
# # 
# depth_factors = ulirschCounts %>%
#   select(matches('NA')) %>%
#   summarise_all(.funs = funs(sum(.) / 1e6)) %>%
#   gather(sample, depth_factor)

# allele = mpra_data$count_data[[1]] %>%
#   pull(allele) %>%
#   {. != 'ref'} %>%
#   as.integer() %>%
#   {. + 1}
# 
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
# 
library(parallel)
# nb_param_estimates = mpra_data %>%
#   unnest %>%
#   gather(sample, count, matches('NA')) %>%
#   group_by(snp_id, allele, sample) %>%
#   nest %>%
#   mutate(count_mean = map_dbl(data, ~mean(.x$count)),
#          nb_fit = mclapply(data, fit_nb, mc.cores = 5),
#          converged = map_lgl(nb_fit, ~.x$convergence == 0)) %>%
#   filter(converged) %>%
#   left_join(sample_depths, by = 'sample') %>%
#   mutate(depth_adj_mean = 1e6 * count_mean / depth,
#          depth_adj_mu_est = map2_dbl(nb_fit, depth, ~1e6 * .x$par[1] / .y),
#          phi_est = map_dbl(nb_fit, ~.x$par[2])) %>%
#   filter(depth_adj_mean > 10)
# 
# nb_param_estimates %>%
#   ggplot(aes(depth_adj_mu_est)) +
#   geom_density(aes(color = sample)) +
#   scale_x_log10()
# 
# library(stringr)
# marg_prior = nb_param_estimates %>%
#   mutate(acid_type = factor(str_extract(sample, 'DNA|RNA'))) %>%
#   group_by(allele, acid_type) %>%
#   summarise(phi_gamma_prior = list(fit_gamma(phi_est)),
#             mu_gamma_prior = list(fit_gamma(depth_adj_mu_est))) %>%
#   ungroup %>%
#   gather(prior_type, gamma_fit, matches('gamma')) %>%
#   mutate(alpha_est = map_dbl(gamma_fit, ~.x$par[1]),
#          beta_est = map_dbl(gamma_fit, ~.x$par[2])) %>%
#   arrange(desc(allele)) # This puts reference alleles first. This is bad practice
# 
# 
# load("/mnt/bigData2/andrew/analysis_data/testing_dat.RData")
# 
# sample_depths = mpra_data %>%
#   unnest %>%
#   select(-snp_id, -allele) %>%
#   gather(sample, count) %>%
#   group_by(sample) %>%
#   summarise(depth = sum(count))
# 
# snp_dat = mpra_data$count_data[[1]] %>%
#   mutate(bc_id = 1:n())
# 
# data_list = list(n_rna_samples = mpra_data$count_data[[1]] %>% select(matches('RNA')) %>% ncol,
#                  n_dna_samples = mpra_data$count_data[[1]] %>% select(matches('DNA')) %>% ncol,
#                  n_barcodes = mpra_data$count_data[[1]] %>% nrow,
#                  rna_counts = mpra_data$count_data[[1]] %>% select(matches('RNA')) %>% as.matrix,
#                  dna_counts = mpra_data$count_data[[1]] %>% select(matches('DNA')) %>% as.matrix,
#                  allele = allele,
#                  rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
#                  dna_depths = depth_factors %>% filter(grepl('DNA', sample)) %>% pull(depth_factor),
#                  #dna_norm_factors = inputs$DNA_norm_factor,
#                  rna_norm_factors = inputs$DNA_norm_factor,
#                  rna_m_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
#                  rna_m_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(beta_est),
#                  rna_p_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(alpha_est), # horrible non-alignment :(
#                  rna_p_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(beta_est),
#                  dna_m_a = marg_prior %>% filter(acid_type == 'DNA', prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
#                  dna_m_b = marg_prior %>% filter(acid_type == 'DNA', prior_type == 'mu_gamma_prior') %>% pull(beta_est),
#                  dna_p_a = marg_prior %>% filter(acid_type == 'DNA', prior_type == 'phi_gamma_prior') %>% pull(alpha_est),
#                  dna_p_b = marg_prior %>% filter(acid_type == 'DNA', prior_type == 'phi_gamma_prior') %>% pull(beta_est))
# 

#### model string ----

# lol this is wrong
# bc_effect_model = '
# data {
#   int<lower=0> n_rna_samples;
#   int<lower=0> n_barcodes; // number for this allele
#   int<lower=0> rna_counts[n_barcodes, n_rna_samples]; 
#   int<lower=1, upper = 2> allele[n_barcodes]; // allele indicator; 1 = ref, 2 = alt
#   real<lower=0> rna_depths[n_rna_samples]; 
#   real<lower=0> rna_norm_factors[n_barcodes];
#   real<lower=0> rna_m_a[2];
#   real<lower=0> rna_m_b[2];
#   real<lower=0> rna_p_a[2];
#   real<lower=0> rna_p_b[2];
# } 
# parameters {
#   vector<lower=0>[2] r_m_i;
#   vector<lower=0>[2] r_p_i;
# }
# model {
# 
#   // with density estimation, alleles would have different priors
#   r_m_i[allele] ~ gamma(rna_m_a[allele], rna_m_b[allele]); // priors on negative binomial parameters
#   r_p_i[allele] ~ gamma(rna_p_a[allele], rna_p_b[allele]); // here, both alleles come from the same prior
# 
#   for (s in 1:n_rna_samples) {
#     for (t in 1:n_barcodes) {
#       rna_counts[allele, s][t] ~ neg_binomial_2(r_m_i[allele] * rna_depths[s] * rna_norm_factors[t], r_p_i[allele]);
#     }
#     
#   }
# 
# }
# generated quantities {
#   real transcription_shift;
#   transcription_shift = log(r_m_i[2]) - log(r_m_i[1]);
# }
# '
# 
# # Divide out average DNA count in TS? TODO
# # Rename to abundance effect
# 
# bc_object = stan_model(model_code = bc_effect_model)

#### test ----

# samp_test = sampling(bc_object,
#                      data = data_list,
#                      chains = 10,
#                      iter = 1300,
#                      warmup = 300,
#                      cores = 10) # ~ 62 seconds


library(coda)
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

run_samp_test = function(count_data, snp_id,
                         save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/',
                         depth_factors,
                         n_cores = 10,
                         tot_samp = 1e4){
  
  snp_dat = count_data %>% 
    mutate(bc_id = 1:n())
  
  # For one barcode take the depth normalized values in each replicate .
  
  dnv = snp_dat %>% 
    select(allele, bc_id, matches('NA')) %>% 
    gather(sample, count, -allele, -bc_id) %>% 
    left_join(depth_factors, by = 'sample') %>% 
    mutate(depth_norm_count = count / depth_factor)
  
  well_represented = dnv %>% 
    filter(grepl('DNA', sample)) %>% 
    group_by(allele, bc_id) %>% 
    summarise(mean_depth_norm = mean(depth_norm_count)) %>% 
    ungroup %>% 
    filter(mean_depth_norm > 10)
  
  wr_counts = well_represented %>% 
    count(allele)
  
  if (any(wr_counts$n < 2) | nrow(well_represented) == 0) {
    return(NA)
  }
  
  dnv %<>% filter(bc_id %in% well_represented$bc_id)
  
  # 2) take the median of these values
  
  medians = dnv %>% 
    mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>% 
    group_by(bc_id, samp_type) %>%
    summarise(med_dnv = median(depth_norm_count)) %>% 
    ungroup %>% 
    left_join(unique(select(dnv, bc_id, allele)),
              by = 'bc_id') 
  
  # 3) apply the geometric mean based procedure above 
  #bc_medians = 
  
  # https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
  gm_mean = function(x, na.rm=TRUE){ 
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
  geo_means = dnv %>% 
    mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>% 
    group_by(samp_type, allele) %>% 
    summarise(geo_mean = gm_mean(depth_norm_count))
  
  bc_norm_factors = medians %>% 
    left_join(geo_means, by = c('samp_type', 'allele')) %>% 
    mutate(bc_norm_factor = (med_dnv / geo_mean)) %>% 
    select(bc_id, samp_type, bc_norm_factor) %>% 
    mutate(samp_type = paste0(samp_type, '_norm_factor')) %>% 
    spread(samp_type, bc_norm_factor)
  
  inputs = snp_dat %>% 
    filter(bc_id %in% well_represented$bc_id) %>% 
    left_join(bc_norm_factors,
              by = 'bc_id') %>% 
    select(allele, bc_id, DNA_norm_factor, RNA_norm_factor, everything())
  
  count_data = snp_dat %>% 
    filter(bc_id %in% well_represented$bc_id)
  
  data_list = list(n_rna_samples = count_data %>% select(matches('RNA')) %>% ncol,
                   n_barcodes = count_data %>% nrow,
                   rna_counts = count_data %>% select(matches('RNA')) %>% as.matrix,
                   allele = count_data %>% mutate(allele_ind = case_when(allele == 'ref' ~ 1, allele == 'mut' ~ 2)) %>% pull(allele_ind), 
                   rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
                   rna_norm_factors = inputs$DNA_norm_factor,
                   rna_m_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
                   rna_m_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(beta_est),
                   rna_p_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(alpha_est), # horrible non-alignment :(
                   rna_p_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(beta_est))
  
  n_samp_per_core = tot_samp / n_cores 
  n_iter = n_samp_per_core + 300
  
  samp_test = sampling(bc_object,
                       data = data_list,
                       chains = n_cores, 
                       iter = n_iter, 
                       warmup = 300,
                       cores = n_cores)
  save(samp_test, data_list,
       file = paste0(save_dir, snp_id %>% gsub(' ', '_', .) %>% gsub('\\/', '-', .), '.RData'))
  
  samp_test %>% 
    rstan::extract() %>% 
    .[['transcription_shift']] %>% 
    mcmc %>% 
    my_HPD
  
}

mean_norm_factor_test = function(count_data, snp_id,
                                 save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/mean_norm_factor/',
                                 depth_factors,
                                 n_cores = 10,
                                 tot_samp = 1e4){
  
  snp_dat = count_data %>% 
    mutate(bc_id = 1:n())
  
  # For one barcode take the depth normalized values in each replicate .
  
  dnv = snp_dat %>% 
    select(allele, bc_id, matches('NA')) %>% 
    gather(sample, count, -allele, -bc_id) %>% 
    left_join(depth_factors, by = 'sample') %>% 
    mutate(depth_norm_count = count / depth_factor)
  
  well_represented = dnv %>% 
    filter(grepl('DNA', sample)) %>% 
    group_by(allele, bc_id) %>% 
    summarise(mean_depth_norm = mean(depth_norm_count)) %>% 
    ungroup %>% 
    filter(mean_depth_norm > 10)
  
  wr_counts = well_represented %>% 
    count(allele)
  
  if (any(wr_counts$n < 2) | nrow(well_represented) == 0) {
    return(NA)
  }
  
  dnv %<>% filter(bc_id %in% well_represented$bc_id)
  
  bc_mean_factors = dnv %>% 
    mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>% 
    group_by(bc_id, samp_type) %>%
    summarise(mean_dnv = mean(depth_norm_count)) %>% # mean depth normalized count by barcode
    ungroup %>% 
    left_join(unique(select(dnv, bc_id, allele)), # attach on allele
              by = 'bc_id') 
  
  samp_means = dnv %>% 
    mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>% 
    group_by(samp_type, allele) %>% 
    summarise(samp_mean = mean(depth_norm_count)) # mean depth normalized count by sample
  
  bc_norm_factors = bc_mean_factors %>% 
    left_join(samp_means, by = c('samp_type', 'allele')) %>% 
    mutate(bc_norm_factor = (mean_dnv / samp_mean)) %>% 
    select(bc_id, samp_type, bc_norm_factor) %>% 
    mutate(samp_type = paste0(samp_type, '_norm_factor')) %>% 
    spread(samp_type, bc_norm_factor)
  
  inputs = snp_dat %>% 
    filter(bc_id %in% well_represented$bc_id) %>% 
    left_join(bc_norm_factors,
              by = 'bc_id') %>% 
    select(allele, bc_id, DNA_norm_factor, RNA_norm_factor, everything())
  
  count_data = snp_dat %>% 
    filter(bc_id %in% well_represented$bc_id)
  
  data_list = list(n_rna_samples = count_data %>% select(matches('RNA')) %>% ncol,
                   n_barcodes = count_data %>% nrow,
                   rna_counts = count_data %>% select(matches('RNA')) %>% as.matrix,
                   allele = count_data %>% mutate(allele_ind = case_when(allele == 'ref' ~ 1, allele == 'mut' ~ 2)) %>% pull(allele_ind), 
                   rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
                   rna_norm_factors = inputs$DNA_norm_factor,
                   rna_m_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
                   rna_m_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(beta_est),
                   rna_p_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(alpha_est), # horrible non-alignment :(
                   rna_p_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(beta_est))
  
  n_samp_per_core = tot_samp / n_cores 
  n_iter = n_samp_per_core + 300
  
  samp_test = sampling(bc_object,
                       data = data_list,
                       chains = n_cores, 
                       iter = n_iter, 
                       warmup = 300,
                       cores = n_cores)
  save(samp_test, data_list,
       file = paste0(save_dir, snp_id %>% gsub(' ', '_', .) %>% gsub('\\/', '-', .), '.RData'))
  
  samp_test %>% 
    rstan::extract() %>% 
    .[['transcription_shift']] %>% 
    mcmc %>% 
    my_HPD
  
}

# bc_effect_tests = mpra_data %>% 
#   mutate(ts_HDI = map2(count_data, snp_id, run_samp_test, depth_factors = depth_factors))
# 
# save(bc_effect_tests, 
#      file = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests.RData')


#### Apply to CD36 results ----
load("~/plateletMPRA/outputs/pcr_validation_pilot/controls_with_counts.RData")
load("~/plateletMPRA/outputs/pcr_validation_pilot/eqtls_with_counts.RData")


controls %<>% mutate(rs = rep(c(paste0('PKRRE', 1:5), paste0('ALAS2', 1:3), 'URUOS', 'HBG2'), each = 80))
eqtls %<>% mutate(rs = map2_chr(rs, mut, ~ifelse(.x == 'rs17154155' & .y == 'T', 'rs17154155_ALT', .x))) # this snp had two alternate alleles

pltMPRA = rbind(eqtls, controls)
pltMPRA %<>% set_colnames(gsub('_L001_R1_001|seqOnly_', 
                               '',
                               names(pltMPRA))) %>% 
  select_all(~gsub('cDNA', 'RNA', gsub('Plasmid', 'DNA', .)))

library(stringr)
pltMPRA %<>%
  select(-seq, allele = type, snp = rs)

cd36MPRA_depth_factors = pltMPRA %>% 
  select(matches('NA')) %>%
  summarise_all(.funs = funs(sum(.) / 1e6)) %>% 
  gather(sample, depth_factor)

##

## estimate priors----

library(parallel)
sample_depths = pltMPRA %>% 
  select(matches('NA')) %>% 
  gather(sample, count) %>% 
  group_by(sample) %>% 
  summarise(depth = sum(count))

nb_param_estimates = pltMPRA %>%
  unnest %>%
  gather(sample, count, matches('[DR]NA')) %>%
  group_by(snp, allele, sample) %>%
  nest %>%
  mutate(count_mean = map_dbl(data, ~mean(.x$count)),
         nb_fit = mclapply(data, fit_nb, mc.cores = 5),
         converged = map_lgl(nb_fit, ~.x$convergence == 0)) %>%
  filter(converged) %>%
  left_join(sample_depths, by = 'sample') %>%
  mutate(depth_adj_mean = 1e6 * count_mean / depth,
         depth_adj_mu_est = map2_dbl(nb_fit, depth, ~1e6 * .x$par[1] / .y),
         phi_est = map_dbl(nb_fit, ~.x$par[2])) %>%
  filter(depth_adj_mean > 10)

nb_param_estimates %>%
  ggplot(aes(depth_adj_mu_est)) +
  geom_density(aes(color = sample)) +
  scale_x_log10()

library(stringr)
marg_prior = nb_param_estimates %>%
  mutate(acid_type = factor(str_extract(sample, 'DNA|RNA'))) %>%
  group_by(allele, acid_type) %>%
  summarise(phi_gamma_prior = list(fit_gamma(phi_est)),
            mu_gamma_prior = list(fit_gamma(depth_adj_mu_est))) %>%
  ungroup %>%
  gather(prior_type, gamma_fit, matches('gamma')) %>%
  mutate(alpha_est = map_dbl(gamma_fit, ~.x$par[1]),
         beta_est = map_dbl(gamma_fit, ~.x$par[2])) %>%
  arrange(desc(allele)) 

## run test ----

# cd36_bc_effect_test = pltMPRA %>% 
#   mutate(allele = tolower(allele)) %>% 
#   group_by(snp) %>% 
#   nest %>% 
#   .[c(88),] %>% 
#   mutate(ts_HDI = map2(data, snp, 
#                        run_samp_test, 
#                        save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/',
#                        depth_factors = cd36MPRA_depth_factors,
#                        n_cores = 10))
# 
# cd36_bc_effect_test = pltMPRA %>% 
#   mutate(allele = tolower(allele)) %>% 
#   group_by(snp) %>% 
#   nest %>% 
#   .[c(82, 88, 90, 43, 8),] %>% 
#   mutate(ts_HDI = map2(data, snp, 
#                        run_samp_test, 
#                        save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/',
#                        depth_factors = cd36MPRA_depth_factors,
#                        n_cores = 10))

get_post_mean = function(snp_id, 
                         dir){
  load(paste0(dir, snp_id, '.RData'))
  
  samp_test %>% 
    rstan::extract() %>% 
    .[['transcription_shift']] %>% 
    mean
}

# cd36_bc_effect_test %<>%
#   mutate(post_mean_ts = map_dbl(snp,
#                                 get_post_mean),
#          functional = map_lgl(ts_HDI, 
#                               ~!between(0, .x[1], .x[2])))


# cd36_bc_effect_test = pltMPRA %>% 
#   mutate(allele = tolower(allele)) %>% 
#   group_by(snp) %>% 
#   nest %>% 
#   mutate(ts_HDI = map2(data, snp, 
#                        run_samp_test, 
#                        save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/',
#                        depth_factors = cd36MPRA_depth_factors,
#                        n_cores = 20)) %>% 
#   mutate(post_mean_ts = map_dbl(snp,
#                                 get_post_mean),
#          functional = map_lgl(ts_HDI, 
#                               ~!between(0, .x[1], .x[2])))
# 
# save(cd36_bc_effect_test, 
#      file = '~/bayesianMPRA/analysis_outputs/cd36_bc_effect_test.RData')
# 
# cd36_mean_bc_effect_test = pltMPRA %>% 
#   mutate(allele = tolower(allele)) %>% 
#   group_by(snp) %>% 
#   nest %>% 
#   mutate(ts_HDI = map2(data, snp, 
#                        mean_norm_factor_test, 
#                        save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/mean_norm_factor/',
#                        depth_factors = cd36MPRA_depth_factors,
#                        n_cores = 18)) %>% 
#   mutate(post_mean_ts = map_dbl(snp,
#                                 get_post_mean),
#          functional = map_lgl(ts_HDI, 
#                               ~!between(0, .x[1], .x[2])))

save(cd36_mean_bc_effect_test, 
     file = '~/bayesianMPRA/analysis_outputs/cd36_mean_bc_effect_test.RData')

make_ts_plot = function(snp_id, 
                        dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/') {
  load(paste0(dir, snp_id, '.RData'))
  
  samp_test %>% 
    rstan::extract() %>% 
    .[['transcription_shift']] %>% 
    data_frame(transcription_shift = .) %>% 
    ggplot(aes(transcription_shift)) + 
    geom_histogram(aes(y = ..density..),
                   bins = 40) + 
    geom_density() + 
    labs(title = paste0(snp_id, ' Barcode effect model TS'),
         subtitle = 'TS = log(alt RNA mean) - log(ref RNA mean) after accounting for depth and barcode')
}
# 
# tmp = mclapply(cd36_mean_bc_effect_test$snp,
#     make_ts_plot, mc.cores = 10)


#### Double normalization model ----

bc_effect_model = '
data {
  int<lower=0> n_rna_samples;
  int<lower=1> n_ref;
  int<lower=1> n_mut;
  int<lower=0> ref_counts[n_ref, n_rna_samples];
  int<lower=0> mut_counts[n_mut, n_rna_samples];
  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> ref_rna_norm_factors[n_ref];
  real<lower=0> mut_rna_norm_factors[n_mut];
  real<lower=0> rna_m_a[2];
  real<lower=0> rna_m_b[2];
  real<lower=0> rna_p_a[2];
  real<lower=0> rna_p_b[2];
} 
parameters {
  vector<lower=0>[2] r_m_i;
  vector<lower=0>[2] r_p_i;
}
model {

  // with density estimation, alleles would have different priors
  for (allele in 1:2) {
    r_m_i[allele] ~ gamma(rna_m_a[allele], rna_m_b[allele]); // priors on negative binomial parameters
    r_p_i[allele] ~ gamma(rna_p_a[allele], rna_p_b[allele]); // here, both alleles come from the same prior
  }
  
  for (s in 1:n_rna_samples) {
    for (t in 1:n_ref) {
      ref_counts[t, s] ~ neg_binomial_2(r_m_i[1] * rna_depths[s] * ref_rna_norm_factors[t], r_p_i[1]);
    }

    for (t in 1:n_mut) {
      mut_counts[t, s] ~ neg_binomial_2(r_m_i[2] * rna_depths[s] * mut_rna_norm_factors[t], r_p_i[2]);
    }
  }

}
generated quantities {
  real transcription_shift;
  transcription_shift = log(r_m_i[2]) - log(r_m_i[1]);
}
'

bc_object = stan_model(model_code = bc_effect_model)


bc_norm_factor_test = function(count_data, snp_id,
                                 save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/mean_norm_factor/',
                                 depth_factors,
                                 n_cores = 10,
                                 tot_samp = 1e4){
  
  snp_dat = count_data %>% 
    mutate(bc_id = 1:n())
  
  # For one barcode take the depth normalized values in each replicate .
  
  dnv = snp_dat %>% # depth_normalized_values
    select(allele, bc_id, matches('[DR]NA')) %>% 
    gather(sample, count, -allele, -bc_id) %>% 
    left_join(depth_factors, by = 'sample') %>% 
    mutate(depth_norm_count = count / depth_factor)
  
  well_represented = dnv %>% 
    filter(grepl('DNA', sample)) %>% 
    group_by(allele, bc_id) %>% 
    summarise(mean_depth_norm = mean(depth_norm_count)) %>% 
    ungroup %>% 
    filter(mean_depth_norm > 10)
  
  wr_counts = well_represented %>% 
    count(allele)
  
  if (any(wr_counts$n < 2) | nrow(well_represented) == 0) {
    return(NA)
  }
  
  dnv %<>% filter(bc_id %in% well_represented$bc_id)
  
  bc_mean_factors = dnv %>% 
    mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>% 
    group_by(bc_id, samp_type) %>%
    summarise(mean_dnv = mean(depth_norm_count)) %>% # mean depth normalized count by barcode
    ungroup %>% 
    left_join(unique(select(dnv, bc_id, allele)), # attach on allele
              by = 'bc_id') 
  
  samp_means = dnv %>% 
    mutate(samp_type = if_else(grepl('DNA', sample), 'DNA', 'RNA')) %>% 
    group_by(samp_type, allele) %>% 
    summarise(samp_mean = mean(depth_norm_count)) %>%  # mean depth normalized count by sample
    ungroup
  
  bc_norm_factors = bc_mean_factors %>% 
    left_join(samp_means, by = c('samp_type', 'allele')) %>% 
    mutate(bc_norm_factor = (mean_dnv / samp_mean)) %>% 
    select(bc_id, samp_type, bc_norm_factor) %>% 
    mutate(samp_type = paste0(samp_type, '_norm_factor')) %>% 
    spread(samp_type, bc_norm_factor)
  
  inputs = snp_dat %>% 
    filter(bc_id %in% well_represented$bc_id) %>% 
    left_join(bc_norm_factors,
              by = 'bc_id') %>% 
    select(allele, bc_id, DNA_norm_factor, RNA_norm_factor, everything())
  
  count_data = snp_dat %>% 
    filter(bc_id %in% well_represented$bc_id)
  
  data_list = list(n_rna_samples = count_data %>% select(matches('RNA')) %>% ncol,
                   n_barcodes = inputs %>% nrow,
                   ref_counts = inputs %>% filter(allele == 'ref') %>%  select(matches('RNA')) %>% select(-matches('norm')) %>% as.matrix,
                   mut_counts = inputs %>% filter(allele == 'mut') %>%  select(matches('RNA')) %>% select(-matches('norm')) %>% as.matrix,
                   n_ref = inputs$allele %>% table() %>% .['ref'],
                   n_mut = inputs$allele %>% table() %>% .['mut'],
                   rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
                   ref_rna_norm_factors = inputs %>% filter(allele == 'ref') %>% pull(DNA_norm_factor),
                   mut_rna_norm_factors = inputs %>% filter(allele == 'mut') %>% pull(DNA_norm_factor),
                   rna_m_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
                   rna_m_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'mu_gamma_prior') %>% pull(beta_est),
                   rna_p_a = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(alpha_est), # horrible non-alignment :(
                   rna_p_b = marg_prior %>% filter(acid_type == 'RNA', prior_type == 'phi_gamma_prior') %>% pull(beta_est))
  
  n_samp_per_core = tot_samp / n_cores 
  n_iter = n_samp_per_core + 300
  
  samp_test = sampling(bc_object,
                       data = data_list,
                       chains = n_cores, 
                       iter = n_iter, 
                       warmup = 300,
                       cores = n_cores)
  save(samp_test, data_list,
       file = paste0(save_dir, snp_id %>% gsub(' ', '_', .) %>% gsub('\\/', '-', .), '.RData'))
  
  samp_test %>% 
    rstan::extract() %>% 
    .[['transcription_shift']] %>% 
    mcmc %>% 
    my_HPD
}

cd36_bc_effect_test = pltMPRA %>% 
  mutate(allele = tolower(allele)) %>% 
  group_by(snp) %>% 
  nest %>% 
  mutate(ts_HDI = map2(data, snp, 
                       bc_norm_factor_test, 
                       save_dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/bc_norm_factor/',
                       depth_factors = cd36MPRA_depth_factors,
                       n_cores = 18)) %>% 
  mutate(post_mean_ts = map_dbl(snp,
                                get_post_mean),
         functional = map_lgl(ts_HDI, 
                              ~!between(0, .x[1], .x[2])))

load_and_get_ts_hdi = function(snp_id){
  load(paste0('/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/bc_norm_factor/', snp_id, '.RData'))
  samp_test %>% 
    rstan::extract() %>% 
    .[['transcription_shift']] %>% 
    mcmc %>% 
    my_HPD
}

cd36_bc_effect_test %<>%
  mutate(ts_HDI = map(snp,
                      load_and_get_ts_hdi),
         functional = map_lgl(ts_HDI, 
                              ~!between(0, .x[1], .x[2])),
         post_mean_ts = map_dbl(snp,
                                get_post_mean,
                                dir = '/mnt/labhome/andrew/bayesianMPRA/analysis_outputs/bc_effect_tests/cd36/bc_norm_factor/'))



save(cd36_mean_bc_effect_test, 
     file = '~/bayesianMPRA/analysis_outputs/cd36_mean_bc_effect_test.RData')





