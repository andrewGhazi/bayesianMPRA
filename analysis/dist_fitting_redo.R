library(tidyverse)
library(parallel)
library(magrittr)

# The plain fitting function fails for all 0's or <2 observations
# load("~/bayesianMPRA/analysis_data/testing_dat.RData")
# mpra_data$count_data[[9]] %>% 
#   gather(sample, count, -allele) %>% 
#   group_by(sample, allele) %>% 
#   summarise(nb_fit = list(fitdist(count, distr = 'nbinom')))

#### Functions ----
fit_nb = function(sample_allele_counts){
  if (sum(sample_allele_counts) == 0) {
    return(NA)
  } else if (length(sample_allele_counts) < 2) {
    return(NA)
  }
  
  fitdistrplus::fitdist(data = sample_allele_counts,
          distr = 'nbinom')
}

est_sample_nb = function(count_dat){
  # count_dat - a tibble with a allele (ref/mut) column and columns of observed MPRA counts in given transfections
  
  count_dat %>% 
    gather(block, count, -allele, -barcode) %>% # Need to add some check or conditional on whether or not it includes a barcode column, probably needs rlang
    group_by(allele, block) %>% 
    summarise(mle_neg_bin = list(fit_nb(count))) %>% 
    ungroup %>% 
    filter(map_lgl(mle_neg_bin, ~class(.x) == 'fitdist')) %>% 
    mutate(nb_mu_est = map_dbl(mle_neg_bin, ~.x$estimate['mu']),
           nb_size_est = map_dbl(mle_neg_bin, ~.x$estimate['size'])) %>% 
    dplyr::select(-mle_neg_bin)
}

filter_poorly_represented = function(count_dat){
  depth_adj_dna = count_dat %>% 
    dplyr::select(barcode, contains('DNA')) %>% 
    gather(sample, count, -barcode) %>% # Need to add some check or conditional on whether or not it includes a barcode column, probably needs rlang
    left_join(total_counts, by = 'sample') %>% 
    mutate(depth_adj_count = 1e6 * count / total_count) %>% 
    group_by(barcode) %>% 
    summarise(dna_depth_adj_mean = mean(depth_adj_count))
  
  good_dna = depth_adj_dna %>% 
    filter(dna_depth_adj_mean > 10**.33)
  
  count_dat %>% 
    filter(barcode %in% good_dna$barcode)
}

### Application to tewhey data -----
load('~/bayesianMPRA/analysis_data/tewhey_subset.RData')

total_counts = tewhey_subset %>% 
  mutate(sample = gsub('ctrl', 'DNA', sample) %>% gsub('HepG2', 'HepG2_RNA', .)) %>% 
  filter(grepl('DNA', sample)) %>% 
  group_by(sample) %>%
  summarise(total_count = sum(count)) 

tewhey_subset %>% 
  mutate(sample = gsub('ctrl', 'DNA', sample) %>% gsub('HepG2', 'HepG2_RNA', .)) %>% 
  left_join(total_counts, by = 'sample') %>%
  mutate(depth_adj_count = count / total_count * 1e6) %>%
  filter(grepl('DNA', sample)) %>% 
  ggplot(aes(depth_adj_count)) +
  geom_density(aes(color = sample)) +
  scale_x_log10() +
  geom_vline(xintercept = 10**.33,
             lty = 2,
             color = 'grey30')

t_sub = tewhey_subset %>% 
  mutate(sample = gsub('ctrl', 'DNA', sample) %>% gsub('HepG2', 'HepG2_RNA', .)) %>% 
  spread(sample, count, fill = 0) %>% 
  group_by(snp_id) %>% 
  nest



t_sub %<>% mutate(data = mclapply(data, filter_poorly_represented, mc.cores = 20),
                  nb_params = mclapply(data, est_sample_nb, mc.cores = 20))

# Assign weights

# Fit Gammas

# Run sampler
