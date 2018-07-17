# Do I need to adjust phi for depth?

library(tidyverse)

load("/mnt/bigData2/andrew/analysis_data/testing_dat.RData")

fit_nb = function(counts){
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

sample_depths = mpra_data %>%
  unnest %>%
  select(-snp_id, -allele) %>%
  gather(sample, count) %>%
  group_by(sample) %>% 
  summarise(depth = sum(count))

nb_fits = mpra_data %>% 
  unnest %>% 
  gather(sample, counts, -snp_id, -allele) %>% 
  group_by(snp_id, allele, sample) %>% 
  summarise(nb_fit = list(fit_nb(counts))) %>% 
  ungroup %>% 
  mutate(phi_estimate = map_dbl(nb_fit, ~.x$par[2])) %>% 
  mutate(converged = map_lgl(nb_fit, ~.x$convergence == 0))

nb_fits %>% 
  filter(converged) %>% 
  left_join(sample_depths) %>% 
  ggplot(aes(depth, phi_estimate)) +
  geom_violin() 

nb_fits %>% 
  filter(converged) %>% filter(grepl('RNA', sample)) %>% 
  left_join(sample_depths) %>% filter(phi_estimate > .01) %>% 
  ggplot(aes(depth, phi_estimate)) +
  geom_jitter(height = 0, alpha = .05) + scale_y_log10() + labs(title = 'Phi estimates by depth, Ulirsch 2016 test set') 

nb_fits %>% 
  filter(converged) %>% 
  mutate(mu_estimate = map_dbl(nb_fit, ~.x$par[1])) %>% 
  filter(grepl('RNA', sample)) %>% 
  left_join(sample_depths) %>%
  filter(phi_estimate > .01) %>% 
  ggplot(aes(depth, mu_estimate)) +
  geom_jitter(height = 0, alpha = .05) + 
  scale_y_log10() + geom_smooth(method = 'lm') + labs(title = 'Mean estimates by depth, Ulirsch 2016 test set')
