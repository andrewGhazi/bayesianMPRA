fit_nb = function(sample_allele_counts){
  if (sum(sample_allele_counts) == 0){
    return(NA)
  }
  
  fitdist(data = sample_allele_counts,
          distr = 'nbinom')
}

mpra_data$count_data[[9]] %>% 
  gather(sample, count, -allele) %>% 
  group_by(sample, allele) %>% 
  summarise(nb_fit = list(fitdist(count, distr = 'nbinom')))


est_sample_nb = function(count_dat){
  # count_dat - a tibble with a allele (ref/mut) column and columns of observed MPRA counts in given transfections
  # uses fitdistrplus::fitdist because MASS::fitdistr was cracking wise at me
  # using a modified version of fitdistrplus::mledist because the regular version can't use non-integer weights (and it was cracking wise)
  
  count_dat %>% 
    gather(block, count, -allele) %>% 
    group_by(allele, block) %>% 
    summarise(mle_neg_bin = list(fit_nb(count))) %>% 
    ungroup %>% 
    filter(map_lgl(mle_neg_bin, ~class(.x) == 'fitdist')) %>% 
    mutate(nb_mu_est = map_dbl(mle_neg_bin, ~.x$estimate['mu']),
           nb_size_est = map_dbl(mle_neg_bin, ~.x$estimate['size'])) %>% 
    dplyr::select(-mle_neg_bin)
}


