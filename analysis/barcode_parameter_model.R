# The barcode effect model runs, but having the barcode effect as a point
# estimate is non-ideal. Let's have it be a parameter that varies in the model.

library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(stringr)


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
  real<lower=0> dna_depths[d_rna_samples]; 
  real<lower=0> ref_rna_norm_factors[n_ref];
  real<lower=0> alt_rna_norm_factors[n_alt];

  real<lower=0> dna_m_a; // ref and alt DNA use the same prior, hence length = 1
  real<lower=0> dna_m_b;
  real<lower=0> dna_p_a;
  real<lower=0> dna_p_b;

  real<lower=0> rna_m_a[2];
  real<lower=0> rna_m_b[2];
  real<lower=0> rna_p_a[2];
  real<lower=0> rna_p_b[2];
} 
parameters {
  vector<lower=0>[n_ref] ref_rna_norm_factors;
  vector<lower=0>[n_alt] alt_rna_norm_factors;

  vector<lower=0>[n_ref] dna_m_ref;
  vector<lower=0>[n_alt] dna_m_alt;
  real<lower=0> dna_p;

  vector<lower=0>[2] rna_m; // rna mean
  vector<lower=0>[2] rna_p; // rna phi aka size
}
model {

  for (t in 1:n_ref) {
    dna_m_ref[t] ~ gamma(dna_m_a, dna_m_p);
  }

  for (t in 1:n_alt) {
    dna_m_alt[t] ~ gamma(dna_m_a, dna_m_p);
  }

  dna_p ~ gamma(dna_p_a, dna_p_b);

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
    rna_p[allele] ~ gamma(rna_p_a[allele], rna_p_b[allele]); // here, both alleles come from the same prior
  }
  
  for (s in 1:n_rna_samples) {
    for (t in 1:n_ref) {
      ref_counts[t, s] ~ neg_binomial_2(rna_m[1] * rna_depths[s] * ref_rna_norm_factors[t], rna_p[1]);
    }

    for (t in 1:n_alt) {
      alt_counts[t, s] ~ neg_binomial_2(rna_m[2] * rna_depths[s] * alt_rna_norm_factors[t], rna_p[2]);
    }
  }

}
generated quantities {
  real transcription_shift;
  transcription_shift = log(rna_m[2]) - log(rna_m[1]);
}
'

bc_object = stan_model(model_code = bc_effect_model)


#### fit the prior ---- 
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

nb_param_estimates = pltMPRA %>%
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
