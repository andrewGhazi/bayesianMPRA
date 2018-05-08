library(tidyverse)
library(magrittr)
library(rstan)

load("/mnt/bigData2/andrew/analysis_data/testing_dat.RData")

# 
# nb_reg_model = '
# data {
#   int<lower=0> n_rna_samples;
#   int<lower=0> n_dna_samples;
#   int<lower=0> n_barcodes; // number in the WHOLE assay, not just per allele
#   int<lower=0> rna_counts[n_barcodes, n_rna_samples]; 
#   int<lower=0> dna_counts[n_barcodes, n_dna_samples]; 
#   int<lower=0> n_alleles;
#   int<lower=0, upper = n_alleles> allele[n_barcodes]; //allele indicator
#   real<lower=0> rna_depths[n_rna_samples]; 
#   real<lower=0> dna_depths[n_rna_samples];
# } 
# parameters {
#   real<lower=0> rna_m_a;
#   real<lower=0> rna_m_b;
#   real<lower=0> dna_m_a;
#   real<lower=0> dna_m_b;
#   real<lower=0> rna_p_a;
#   real<lower=0> rna_p_b;
#   real<lower=0> dna_p_a;
#   real<lower=0> dna_p_b;
# 
#   real r_m_i[n_alleles,n_rna_samples];
#   real r_p_i[n_alleles,n_rna_samples];
#   real d_m_i[n_alleles,n_rna_samples];
#   real d_p_i[n_alleles,n_rna_samples];
# }
# model {
# 
#   
#   rna_m_a ~ gamma(20,20);
#   rna_m_b ~ gamma(20,20);
#   rna_p_a ~ gamma(20,20);
#   rna_p_b ~ gamma(20,20);
# 
#   dna_m_a ~ gamma(20,20);
#   dna_m_b ~ gamma(20,20);
#   dna_p_a ~ gamma(20,20);
#   dna_p_b ~ gamma(20,20);
# 
#   for (s in 1:n_rna_samples) {
#     r_m_i[allele,s] ~ gamma(rna_m_a, rna_m_b);
#     r_p_i[allele,s] ~ gamma(rna_p_a, rna_p_b);
#     rna_counts[allele, s] ~ neg_binomial_2(r_m_i[allele], r_p_i[allele]);
#   }
# 
#   for (s in 1:n_dna_samples){
#     d_m_i[allele,s] ~ gamma(dna_m_a, dna_m_b);
#     d_p_i[allele,s] ~ gamma(dna_p_a, dna_p_b);
#     dna_counts[allele, s] ~ neg_binomial_2(d_m_i[allele], d_p_i[allele]);
#   }
# }
# '
# 
# nb_reg_object = stan_model(model_code = nb_reg_model)
# 
# data_list = list()

#### basal mean ----

nb_basal_model = '
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=0> n_barcodes; // number in the WHOLE assay, not just per allele
  int<lower=0> rna_counts[n_barcodes, n_rna_samples]; 
  int<lower=0> dna_counts[n_barcodes, n_dna_samples]; 
  int<lower=0> n_alleles;
  int<lower=0, upper = n_alleles> allele[n_barcodes]; //allele indicator
  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> dna_depths[n_dna_samples];
} 
parameters {
  real<lower=0> rna_m_a;
  real<lower=0> rna_m_b;
  real<lower=0> dna_m_a;
  real<lower=0> dna_m_b;
  real<lower=0> rna_p_a;
  real<lower=0> rna_p_b;
  real<lower=0> dna_p_a;
  real<lower=0> dna_p_b;

  vector<lower=0>[n_alleles] r_m_i;
  vector<lower=0>[n_alleles] r_p_i;
  vector<lower=0>[n_alleles] d_m_i;
  vector<lower=0>[n_alleles] d_p_i;
}
model {
  rna_m_a ~ gamma(1,1); // priors on allele level gamma parameters for RNA and DNA
  rna_m_b ~ gamma(1,1);
  rna_p_a ~ gamma(1,1);
  rna_p_b ~ gamma(1,1);
  
  dna_m_a ~ gamma(1,1);
  dna_m_b ~ gamma(1,1);
  dna_p_a ~ gamma(1,1);
  dna_p_b ~ gamma(1,1);
  
  r_m_i[allele] ~ gamma(rna_m_a, rna_m_b); // priors on negative binomial parameters
  r_p_i[allele] ~ gamma(rna_p_a, rna_p_b);

  d_m_i[allele] ~ gamma(dna_m_a, dna_m_b);
  d_p_i[allele] ~ gamma(dna_p_a, dna_p_b);

  for (s in 1:n_rna_samples) {
    rna_counts[allele, s] ~ neg_binomial_2(r_m_i[allele] * rna_depths[s], r_p_i[allele]);
  }

  for (s in 1:n_dna_samples){
    dna_counts[allele, s] ~ neg_binomial_2(d_m_i[allele] * dna_depths[s], d_p_i[allele]);
  }
}
'

nb_basal_object = stan_model(model_code = nb_basal_model)

load('/mnt/bigData2/andrew/analysis_data/testing_dat_full.RData')
depth_factors = ulirschCounts %>%
  select(matches('NA')) %>%
  summarise_all(.funs = funs(sum(.) / 1e6)) %>% 
  gather(sample, depth_factor)

id_vec = mpra_data %>% 
  unnest() %>% 
  unite(col = allele_id, snp_id, allele) %>% 
  pull(allele_id) %>% 
  unique() %>% 
  set_names(1:(2*910), nm = .)

allele = mpra_data %>% 
  unnest() %>% 
  unite(col = allele_id, snp_id, allele) %>% 
  mutate(allele_num = id_vec[allele_id]) %>% 
  pull(allele_num)

data_list = list(n_rna_samples = mpra_data$count_data %>% bind_rows %>% select(matches('RNA')) %>% ncol,
                 n_dna_samples = mpra_data$count_data %>% bind_rows %>% select(matches('DNA')) %>% ncol,
                 n_barcodes = mpra_data$count_data %>% bind_rows %>% nrow,
                 rna_counts = mpra_data$count_data %>% bind_rows %>% select(matches('RNA')) %>% as.matrix,
                 dna_counts = mpra_data$count_data %>% bind_rows %>% select(matches('DNA')) %>% as.matrix, 
                 n_alleles = nrow(mpra_data) * 2,
                 allele = allele, 
                 rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
                 dna_depths = depth_factors %>% filter(grepl('DNA', sample)) %>% pull(depth_factor))

# real<lower=0> rna_m_a;
# real<lower=0> rna_m_b;
# real<lower=0> dna_m_a;
# real<lower=0> dna_m_b;
# real<lower=0> rna_p_a;
# real<lower=0> rna_p_b;
# real<lower=0> dna_p_a;
# real<lower=0> dna_p_b;
# 
# vector<lower=0>[n_alleles] r_m_i;
# vector<lower=0>[n_alleles] r_p_i;
# vector<lower=0>[n_alleles] d_m_i;
# vector<lower=0>[n_alleles] d_p_i;

mean_inits = mpra_data %>%
  unnest() %>%
  group_by(snp_id, allele) %>% 
  summarise_all(mean) %>%
  gather(sample, mean_count, matches('NA')) %>%
  left_join(depth_factors, by = 'sample') %>%
  mutate(mean_i = mean_count / depth_factor) %>%
  ungroup %>% 
  group_by(snp_id, allele) %>%
  summarise(dna_mean_i = mean(mean_i[grepl('DNA', sample)]),
            rna_mean_i = mean(mean_i[grepl('RNA', sample)])) %>%
  ungroup

ordered_mean_inits = mpra_data %>% 
  unnest() %>% 
  group_by(snp_id, allele) %>% 
  nest %>% 
  left_join(mean_inits, by = c('snp_id', 'allele')) # we do a join to make sure they get passed in the correct order

give_init_list = function(){
  list(rna_m_a = runif(1, min = 0, max = 3), # trying to mimic default stan initialization
       rna_m_b = runif(1, min = 0, max = 3),
       dna_m_a = runif(1, min = 0, max = 3),
       dna_m_b = runif(1, min = 0, max = 3),
       rna_p_a = runif(1, min = 0, max = 3),
       rna_p_b = runif(1, min = 0, max = 3),
       dna_p_a = runif(1, min = 0, max = 3),
       dna_p_b = runif(1, min = 0, max = 3),
       r_m_i = ordered_mean_inits$rna_mean_i, # initialize at observed means
       r_p_i = runif(1820, min = 0, max = 3),
       d_m_i = ordered_mean_inits$dna_mean_i,
       d_p_i = runif(1820, min = 0, max = 3))
}




# vb_test = vb(object = nb_basal_object,
#              data = data_list)
# 
# vb_test %>% rstan::extract() %>% map(as_tibble) %>% map2(., names(.), ~set_names(.x, paste(.y, names(.x), sep = '_')))


test_samp = sampling(nb_basal_object,
                     data = data_list,
                     cores = 4,
                     chains = 4,
                     iter = 2000,
                     warmup = 200,
                     init = give_init_list)

test_samp2 = sampling(nb_basal_object,
                     data = data_list,
                     cores = 4,
                     chains = 4,
                     iter = 2000,
                     warmup = 200,
                     init = give_init_list)
save(test_samp, file = '/mnt/bigData2/andrew/analysis_outputs/test_samp.RData')
save(test_samp2, file = '/mnt/bigData2/andrew/analysis_outputs/test_samp2.RData')


#### basal NON-hierarchical ----
# This actually retains one layer of 

nonhier_basal_model = '
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=0> n_barcodes; // number for this allele
  int<lower=0> rna_counts[n_barcodes, n_rna_samples]; 
  int<lower=0> dna_counts[n_barcodes, n_dna_samples]; 
  int<lower=1, upper = 2> allele[n_barcodes]; // allele indicator
  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> dna_depths[n_dna_samples];
} 
parameters {
  real<lower=0> rna_m_a;
  real<lower=0> rna_m_b;
  real<lower=0> dna_m_a;
  real<lower=0> dna_m_b;

  real<lower=0> rna_p_a;
  real<lower=0> rna_p_b;
  real<lower=0> dna_p_a;
  real<lower=0> dna_p_b;

  vector<lower=0>[2] r_m_i;
  vector<lower=0>[2] r_p_i;
  vector<lower=0>[2] d_m_i;
  vector<lower=0>[2] d_p_i;
}
model {
  rna_m_a ~ gamma(1,1); // priors on allele level gamma parameters for RNA and DNA
  rna_m_b ~ gamma(1,1);
  rna_p_a ~ gamma(1,1);
  rna_p_b ~ gamma(1,1);

  dna_m_a ~ gamma(1,1);
  dna_m_b ~ gamma(1,1);
  dna_p_a ~ gamma(1,1);
  dna_p_b ~ gamma(1,1);

  // with density estimation, alleles would have different priors
  r_m_i[allele] ~ gamma(rna_m_a, rna_m_b); // priors on negative binomial parameters
  r_p_i[allele] ~ gamma(rna_p_a, rna_p_b); // here, both alleles come from the same prior

  d_m_i[allele] ~ gamma(dna_m_a, dna_m_b);
  d_p_i[allele] ~ gamma(dna_p_a, dna_p_b);

  for (s in 1:n_rna_samples) {
    rna_counts[allele, s] ~ neg_binomial_2(r_m_i[allele] * rna_depths[s], r_p_i[allele]);
  }

  for (s in 1:n_dna_samples){
    dna_counts[allele, s] ~ neg_binomial_2(d_m_i[allele] * dna_depths[s], d_p_i[allele]);
  }
}
generated quantities {
  real transcription_shift;
  transcription_shift = log(r_m_i[2] / d_m_i[2]) - log(r_m_i[1] / d_m_i[1]);
}
'


nonhier_basal_object = stan_model(model_code = nonhier_basal_model)

load('/mnt/bigData2/andrew/analysis_data/testing_dat_full.RData')
depth_factors = ulirschCounts %>%
  select(matches('NA')) %>%
  summarise_all(.funs = funs(sum(.) / 1e6)) %>% 
  gather(sample, depth_factor)

allele = mpra_data$count_data[[1]] %>% 
  pull(allele) %>% 
  {. != 'ref'} %>% 
  as.integer() %>% 
  {. + 1}

fit_nb = function(counts){
  fn_to_min = function(param_vec){
    # param_vec[1] nb mean
    # param_vec[2] nb size
    -sum(dnbinom(counts,
                 mu = param_vec[1], 
                 size = param_vec[2], 
                 log = TRUE))
  }
  
  optim(par = c(100, 1),
        fn_to_min,
        lower = c(0,0))
}

mpra_data %>% 
  unnest %>%
  gather(sample, count, matches('NA')) %>% 
  group_by(snp_id, allele, sample) %>% 
  summarise(nb_fit = list(fit_nb(count)))

data_list = list(n_rna_samples = mpra_data$count_data[[1]] %>% select(matches('RNA')) %>% ncol,
                 n_dna_samples = mpra_data$count_data[[1]] %>% select(matches('DNA')) %>% ncol,
                 n_barcodes = mpra_data$count_data[[1]] %>% nrow,
                 rna_counts = mpra_data$count_data[[1]] %>% select(matches('RNA')) %>% as.matrix,
                 dna_counts = mpra_data$count_data[[1]] %>% select(matches('DNA')) %>% as.matrix, 
                 allele = allele, 
                 rna_depths = depth_factors %>% filter(grepl('RNA', sample)) %>% pull(depth_factor),
                 dna_depths = depth_factors %>% filter(grepl('DNA', sample)) %>% pull(depth_factor))

nonhier_test2 = sampling(object = nonhier_basal_object,
         data = data_list,
         cores = 4)
