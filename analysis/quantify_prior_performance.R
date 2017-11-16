# I want to see how well the conditional prior estimates improve upon the 
# marginal prior. To do this, I will take the mean of the gamma hyper prior
# (which is derived from similar variants) and see how well it correlates with
# mean of the estimated negative binomial distribution.

library(tidyverse)
library(magrittr)
library(fitdistrplus)
library(ghzutils)
#source('~/utility/scripts/utilityFunctions.R')

#### Estimate marginal neg bin distribution on counts in each sample ----
load("~/bayesianMPRA/outputs/varInfoWithNegBinAndGammaParams.RData")

marg_priors = varInfo %>% 
  dplyr::select(construct, negBinParams) %>% 
  unnest %>% 
  na.omit %>% 
  group_by(block) %>% 
  summarise(mu_gamma_marg_hyperprior = list(fitdist(data = muEst,
                                            distr = 'gamma'))) %>% 
  mutate(gamma_mean = map_dbl(mu_gamma_marg_hyperprior, 
                              ~.x$estimate['shape']/.x$estimate['rate']))

#### compare to estimated priors ----

eval_prior_performance = function(neg_bin_param_set, RNA_gamma_param_set){
  RNA_gamma_param_set %>% 
    left_join(marg_priors, by = 'block') %>%
    left_join(neg_bin_param_set, by = c('type', 'block')) %>% 
    rename(mu_gamma_conditional_hyperprior = muGammaHyperPriors) %>% 
    mutate(conditional_prior_dens = map2_dbl(muEst, mu_gamma_conditional_hyperprior, 
                                             ~dgamma(.x, 
                                                     shape = .y['shape'],
                                                     rate = .y['rate'])),
           marg_prior_dens = map2_dbl(muEst, mu_gamma_marg_hyperprior, 
                                      ~dgamma(.x, 
                                              shape = .y$estimate['shape'],
                                              rate = .y$estimate['rate'])),
           prior_ratio = conditional_prior_dens / marg_prior_dens) %>% 
    dplyr::select(type, block, prior_ratio)
}

prior_gain = varInfo %>% 
  dplyr::select(construct, transcriptionalShift, negBinParams, RNAgammaParams) %>% 
  mutate(conditional_improvement = map2(negBinParams, RNAgammaParams, eval_prior_performance))
  
prior_gain %>% 
  dplyr::select(construct, transcriptionalShift, conditional_improvement) %>% 
  unnest %>%
  mutate(log_prior_ratio = log(prior_ratio)) %>% 
  na.omit %>%
  pull(log_prior_ratio) %>% 
  summary()

prior_gain %>% dplyr::select(construct, conditional_improvement) %>% unnest %>% mutate(log_prior_ratio = log(prior_ratio)) %>% na.omit %>% ggplot(aes(log_prior_ratio)) + geom_histogram(bins = 40) + geom_rug() + xlim(-2.5,2.5)

prior_gain %>% 
  dplyr::select(construct, transcriptionalShift, conditional_improvement) %>% 
  unnest %>% ggplot(aes(transcriptionalShift, prior_ratio)) + 
  geom_point(alpha = .1) +
  scale_y_log10(limits = c(exp(-2.5),exp(2.5))) + 
  theme_pres() + 
  geom_hline(yintercept = 1, 
             lty = 2,
             color = 'grey60') + 
  labs(y = 'ML Parameter estimate density ratio \nunder conditional and marginal prior')
ggsave('~/bayesianMPRA/analysis_outputs/plots/prior_gain.png')

prior_gain %>% 
  dplyr::select(construct, transcriptionalShift, conditional_improvement) %>% 
  unnest %>% 
  pull(prior_ratio) %>% 
  {. > 1} %>% 
  na.omit %>% 
  {sum(.) / length(.)}
# [1] 0.7050819
# So the conditional prior agrees more with the maximum likelihood estimate 70% of the time

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -66.6400  -0.0456   0.1683   0.1817   0.3764 482.2000 

# varInfo %>%
#   pull1(RNAgammaParams) %>%
#   mutate(conditional_prior_mean = map_dbl(muGammaHyperPriors, ~.x['shape'] / .x['rate'])) %>%
#   left_join(varInfo %>% pull1(negBinParams), by = c('type', 'block')) %>% 
#   ggplot(aes(muEst, conditional_prior_mean)) + 
#   geom_point()

