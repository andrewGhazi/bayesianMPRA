# The barcode effect model runs, but having the barcode effect as a point
# estimate is non-ideal. Let's have it be a parameter that varies in the model.

library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(stringr)
select = dplyr::select


#### prior fitting functions ----
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


#### new stan model ----

bc_effect_model = '
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=1> n_ref; // number of reference barcodes
  int<lower=1> n_alt; // number of alternate barcodes
  int<lower=0> ref_dna_counts[n_ref, n_dna_samples];
  int<lower=0> alt_dna_counts[n_alt, n_dna_samples];
  int<lower=0> ref_rna_counts[n_ref, n_rna_samples];
  int<lower=0> alt_rna_counts[n_alt, n_rna_samples];

  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> dna_depths[n_dna_samples]; 

  real<lower=0> dna_m_a; // ref and alt DNA use the same prior, hence length = 1
  real<lower=0> dna_m_b; // dna_m are basically the barcode effects here

  real<lower=0> dna_p_a; // dna_p priors are estimated differently
  real<lower=0> dna_p_b; // I should rename these variables

  real<lower=0> rna_m_a[2];
  real<lower=0> rna_m_b[2];
  real<lower=0> rna_p_a[2];
  real<lower=0> rna_p_b[2];
} 
parameters {
  vector<lower=0>[n_ref] dna_m_ref;
  vector<lower=0>[n_alt] dna_m_alt;
  real<lower=0> dna_p;

  vector<lower=0>[2] rna_m; // rna mean
  vector<lower=0>[2] rna_p; // rna phi aka size
}
model {

  for (t in 1:n_ref) {
    dna_m_ref[t] ~ gamma(dna_m_a, dna_m_b);
  }

  for (t in 1:n_alt) {
    dna_m_alt[t] ~ gamma(dna_m_a, dna_m_b);
  }

  dna_p ~ gamma(dna_p_a, dna_p_b);

  // so here the individual counts come from a barcode-specific distribution, 
  // so the gamma prior on dna_m_a/b above needs to be fit on the 
  // depth adjusted counts themselves, not the mean
  for (s in 1:n_dna_samples) {
    for (t in 1:n_ref) {
      ref_dna_counts[t,s] ~ neg_binomial_2(dna_m_ref[t] * dna_depths[s], dna_p);
    }
    
    for (t in 1:n_alt) {
      alt_dna_counts[t,s] ~ neg_binomial_2(dna_m_alt[t] * dna_depths[s], dna_p);
    }
  }
  
  for (allele in 1:2) {
    rna_m[allele] ~ gamma(rna_m_a[allele], rna_m_b[allele]); // priors on negative binomial parameters
    rna_p[allele] ~ gamma(rna_p_a[allele], rna_p_b[allele]); // here, alleles have separate priors
  }
  
  for (s in 1:n_rna_samples) {
    for (t in 1:n_ref) {
      ref_rna_counts[t, s] ~ neg_binomial_2(rna_m[1] * rna_depths[s] * dna_m_ref[t], rna_p[1]);
    }

    for (t in 1:n_alt) {
      alt_rna_counts[t, s] ~ neg_binomial_2(rna_m[2] * rna_depths[s] * dna_m_alt[t], rna_p[2]);
    }
  }

}
generated quantities {
  real ref_act;
  real alt_act;
  real transcription_shift;
  
  ref_act = log(rna_m[1]);
  alt_act = log(rna_m[2]);
  transcription_shift = alt_act - ref_act;
}
'

bc_object = stan_model(model_code = bc_effect_model)


#### fit the priors ---- 
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

# cd36MPRA_depth_factors = pltMPRA %>% 
#   select(matches('NA')) %>%
#   summarise_all(.funs = funs(sum(.) / 1e6)) %>% 
#   gather(sample, depth_factor)


sample_depths = pltMPRA %>% 
  select(matches('NA')) %>% 
  gather(sample, count) %>% 
  group_by(sample) %>% 
  summarise(depth = sum(count))

depth_norm_dna_counts = pltMPRA %>%
  gather(sample, count, matches('[DR]NA')) %>%
  filter(grepl('DNA', sample)) %>% 
  left_join(sample_depths, by = 'sample') %>% 
  mutate(depth_norm_count = 1e6 * count / depth) %>% 
  filter(depth_norm_count > 10)

nb_param_estimates = pltMPRA %>%
  filter(totindex %in% depth_norm_dna_counts$totindex) %>% # only well_represented barcodes
  gather(sample, count, matches('[DR]NA')) %>%
  filter(grepl('[DR]NA', sample)) %>% 
  group_by(snp, allele, sample) %>%
  nest %>%
  mutate(count_mean = map_dbl(data, ~mean(.x$count)),
         nb_fit = mclapply(data, fit_nb, mc.cores = 5),
         converged = map_lgl(nb_fit, ~.x$convergence == 0)) %>%
  filter(converged) %>%
  left_join(sample_depths, by = 'sample') %>%
  mutate(depth_adj_mean = 1e6 * count_mean / depth,
         depth_adj_mu_est = map2_dbl(nb_fit, depth, ~1e6 * .x$par[1] / .y),
         phi_est = map_dbl(nb_fit, ~.x$par[2])) 

scaled_dgamma = function(x, shape, rate, scale){
  scale * dgamma(x, shape = shape, rate = rate)
}

nb_param_estimates %>% # the scaled param estimates show roughly the same distribution by sample
  ggplot(aes(depth_adj_mu_est)) +
  geom_density(aes(color = sample)) +
  scale_x_log10() +
  stat_function(fun = scaled_dgamma,
                args = list(shape = 5.00, rate = .0332, scale = 50),
                color = 'red') + 
  annotate(label = 'DNA prior', color = 'red', x = 1000, geom = 'text', y = 1) + 
  annotate(label = 'Ref RNA prior', color = 'blue', x = 1000, geom = 'text', y = .875) + 
  annotate(label = 'Mut RNA prior', color = 'darkorchid1', x = 1000, geom = 'text', y = .75) + 
  stat_function(fun = scaled_dgamma,
                args = list(shape = .585, rate = .00309, scale = 50),
                color = 'blue') + 
  stat_function(fun = scaled_dgamma,
                args = list(shape = 1.06, rate = .0107, scale = 50),
                color = 'darkorchid1') 

nb_param_estimates %>% # the scaled param estimates show roughly the same distribution by sample
  ggplot(aes(phi_est)) +
  geom_density(aes(color = sample)) +
  scale_x_log10() # hmm, the 22 cycle RNA sample seems to have higher dispersions

# depth_norm_dna_counts %>% 
#   ggplot(aes(depth_norm_count)) +
#   geom_histogram(aes(y = ..density..)) +
#   stat_function(fun = dgamma, 
#                 args = list(shape = dna_gamma_estimates$par[1], 
#                             rate = dna_gamma_estimates$par[2]))
# the gamma fits well, stat_function + scale_x_log10 do not play nicely

library(stringr)
marg_rna_prior = nb_param_estimates %>%
  mutate(acid_type = factor(str_extract(sample, 'DNA|RNA'))) %>%
  group_by(allele, acid_type) %>%
  summarise(phi_gamma_prior = list(fit_gamma(phi_est)),
            mu_gamma_prior = list(fit_gamma(depth_adj_mu_est))) %>%
  ungroup %>%
  gather(prior_type, gamma_fit, matches('gamma')) %>%
  mutate(alpha_est = map_dbl(gamma_fit, ~.x$par[1]),
         beta_est = map_dbl(gamma_fit, ~.x$par[2])) %>%
  arrange(desc(allele)) %>% 
  filter(acid_type == 'RNA')


# This is needed for the prior on Phi for the DNA counts
marg_dna_prior = nb_param_estimates %>%
  mutate(acid_type = factor(str_extract(sample, 'DNA|RNA'))) %>%
  group_by(acid_type) %>% # only group by acid_type here, not allele
  summarise(phi_gamma_prior = list(fit_gamma(phi_est)),
            mu_gamma_prior = list(fit_gamma(depth_adj_mu_est))) %>%
  ungroup %>%
  gather(prior_type, gamma_fit, matches('gamma')) %>%
  mutate(alpha_est = map_dbl(gamma_fit, ~.x$par[1]),
         beta_est = map_dbl(gamma_fit, ~.x$par[2])) %>%
  filter(acid_type == 'DNA') 



# components of data{} input: 
#
# int<lower=0> n_rna_samples;
# int<lower=0> n_dna_samples;
# int<lower=1> n_ref; // number of reference barcodes
# int<lower=1> n_alt; // number of alternate barcodes
# int<lower=0> ref_dna_counts[n_ref, n_dna_samples];
# int<lower=0> alt_dna_counts[n_alt, n_dna_samples];
# int<lower=0> ref_rna_counts[n_ref, n_rna_samples];
# int<lower=0> alt_rna_counts[n_alt, n_rna_samples];
# 
# real<lower=0> rna_depths[n_rna_samples]; 
# real<lower=0> dna_depths[n_dna_samples]; 
# 
# real<lower=0> dna_m_a; // ref and alt DNA use the same prior, hence length = 1
# real<lower=0> dna_m_b;
# real<lower=0> dna_p_a;
# real<lower=0> dna_p_b;
# 
# real<lower=0> rna_m_a[2];
# real<lower=0> rna_m_b[2];
# real<lower=0> rna_p_a[2];
# real<lower=0> rna_p_b[2];

# inputs = pltMPRA %>% 
#   filter(totindex %in% depth_norm_dna_counts$totindex) %>% 
#   group_by(snp) %>% 
#   nest

depth_factors = sample_depths %>% 
  mutate(depth_factor = depth / 1e6)


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

depth_factors = sample_depths %>% 
  mutate(depth_factor = depth / 1e6)


#### function to fit the model ----
run_full_model = function(snp_id, 
                          snp_dat,
                          depth_factors,
                          marg_dna_prior,
                          marg_rna_prior, 
                          n_chains,
                          tot_samp,
                          n_cores,
                          n_warmup){
  

  
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
                   dna_m_a = marg_dna_prior %>% filter(prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
                   dna_m_b = marg_dna_prior %>% filter(prior_type == 'mu_gamma_prior') %>% pull(beta_est),
                   dna_p_a = marg_dna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(alpha_est), 
                   dna_p_b = marg_dna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(beta_est),
                   rna_m_a = marg_rna_prior %>% filter(prior_type == 'mu_gamma_prior') %>% pull(alpha_est),
                   rna_m_b = marg_rna_prior %>% filter(prior_type == 'mu_gamma_prior') %>% pull(beta_est),
                   rna_p_a = marg_rna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(alpha_est), 
                   rna_p_b = marg_rna_prior %>% filter(prior_type == 'phi_gamma_prior') %>% pull(beta_est))
  
  n_samp_per_core = tot_samp / n_chains 
  n_iter = n_samp_per_core + n_warmup
  
  sampler_res = sampling(bc_object,
                       data = data_list,
                       chains = n_chains, 
                       iter = n_iter, 
                       warmup = n_warmup,
                       cores = n_cores)
  
  save(sampler_res,
       file = paste0('/mnt/bigData2/andrew/analysis_outputs/full_model_results/', snp_id, '.RData'))
  
  sampler_res %>% 
    rstan::extract(pars = 'transcription_shift') %>% 
    .[['transcription_shift']] %>% 
    coda::mcmc() %>% 
    my_HPD()
}

full_model_res = pltMPRA %>%
  filter(totindex %in% depth_norm_dna_counts$totindex) %>%
  group_by(snp) %>%
  nest(.key = snp_dat) %>%
  mutate(ts_hdi = mcmapply(run_full_model,
                           snp,
                           snp_dat,
                           MoreArgs = list(depth_factors = depth_factors,
                                           marg_dna_prior = marg_dna_prior,
                                           marg_rna_prior = marg_rna_prior,
                                           tot_samp = 1e4,
                                           n_cores = 1,
                                           n_warmup = 500,
                                           n_chains = 4),
                           SIMPLIFY = FALSE,
                           USE.NAMES = FALSE,
                           mc.cores = 19))

# load('/mnt/bigData2/andrew/analysis_outputs/full_model_results/full_model_res.RData')

get_hdi = function(snp_id, prob){
  load(paste0('/mnt/bigData2/andrew/analysis_outputs/full_model_results/', snp_id, '.RData'))

  sampler_res %>%
    rstan::extract(pars = 'transcription_shift') %>%
    .[['transcription_shift']] %>%
    coda::mcmc() %>%
    my_HPD(prob = prob)
}

get_post_mean = function(snp_id){
  load(paste0('/mnt/bigData2/andrew/analysis_outputs/full_model_results/', snp_id, '.RData'))
  sampler_res %>%
    rstan::extract(pars = 'transcription_shift') %>%
    .[['transcription_shift']] %>%
    mean
}

full_model_res %<>%
  mutate(post_mean = unlist(mclapply(snp, get_post_mean, mc.cores = 19)),
         lower = map_dbl(ts_hdi, ~.x[1]),
         upper = map_dbl(ts_hdi, ~.x[2]),
         functional = map_lgl(ts_hdi, ~!between(0, .x[1], .x[2])),
         hdi_99 = mclapply(snp,
                           get_hdi, 
                           prob = .99,
                           mc.cores = 19),
         functional_99 = map_lgl(hdi_99, ~!between(0, .x[1], .x[2])))


save(full_model_res,
     file = '/mnt/bigData2/andrew/analysis_outputs/full_model_results/full_model_res.RData')

# full_model_res %>%
#   filter(functional)
# 
# # compare to activity t-test results ----
# load('~/plateletMPRA/outputs/pcr_validation_pilot/activities.RData')
# activity_tests = activities %>%
#   group_by(snp) %>% # on each SNP
#   summarise(t_test = list(t.test(activity[allele == 'Mut'], activity[allele == 'Ref'])), # run an t-test on activity levels between alleles
#             u_test = list(wilcox.test(activity[allele == 'Mut'], activity[allele == 'Ref']))) %>%  #also a U-test
#   mutate(tidy_t = map(t_test, broom::tidy), # extract estimates and significance levels
#          transcription_shift = map_dbl(tidy_t, ~.x$estimate),
#          t_p = map_dbl(tidy_t, ~.x$p.value),
#          t_q = p.adjust(t_p, method = 'fdr'),
#          tidy_u = map(u_test, broom::tidy),
#          u_p = map_dbl(tidy_u, ~.x$p.value),
#          u_q = p.adjust(u_p, method = 'fdr')) %>%
#   arrange(t_q)
# 
# full_model_res %>%
#   left_join(activity_tests %>% select(snp, transcription_shift, t_q), by = 'snp') %>%
#   mutate(t_hit = t_q < .05) %>%
#   rename(functional_95 = functional) %>%
#   ggplot(aes(transcription_shift, post_mean)) +
#   geom_point(aes(color = t_hit,
#                  shape = functional_95),
#              size = 3,
#              alpha = .5) +
#   geom_abline(slope = 1, intercept = 0,
#               lty = 2, color = 'grey60') +
#   labs(title = 'Difference in allele means vs. Posterior means')
# 
# sampler_res %>%
#   rstan::extract(pars = 'rna_m') %>% .$rna_m %>%  as_tibble %>% set_colnames(c('ref rna_m', 'alt rna_m')) %>% gather(param, val) %>% ggplot(aes(val)) + geom_histogram(bins = 100, color = 'black') + facet_grid(param ~ ., scales = 'free_x')
# 
