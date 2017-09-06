library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(fitdistrplus)
source('~/bayesianMPRA/R/mledistModified.R')

load('~/bayesianMPRA/outputs/objects_for_isolated_run.RData')


sampler_results = varInfo %>% 
  group_by(construct) %>% 
  nest %>% 
  mutate(sampler_result = mclapply(data, run_sampler, mc.cores = 20))

save(sampler_results,
     file = '~/bayesianMPRA/outputs/sampler_results.RData')
