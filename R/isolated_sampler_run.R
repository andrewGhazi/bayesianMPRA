library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(fitdistrplus)
library(preprocessCore)
library(coda)
source('~/bayesianMPRA/R/mledistModified.R')

load('~/bayesianMPRA/outputs/objects_for_isolated_run.RData')


# sampler_results = varInfo %>% 
#   group_by(construct) %>% 
#   nest %>% 
#   mutate(sampler_result = mclapply(data, run_sampler, mc.cores = 20))
# 
# save(sampler_results,
#      file = '~/bayesianMPRA/outputs/sampler_results.RData')

constructResults = sampler_results$sampler_result

varInfo = sampler_results %>% 
  dplyr::select(-sampler_result) %>% 
  unnest %>% 
  mutate(sampler_result = constructResults)

quant_norm_parameter = function(param_name, sample_df){
  sample_df %>% 
    dplyr::select(contains(param_name)) %>% 
    as.matrix %>% 
    normalize.quantiles %>% 
    rowMeans %>% # Still not sure that taking the mean across samples is the right thing to do
    as.tibble %>% 
    set_names(param_name)
}

get_snp_TS_samples = function(sampler_res){
  sample_df = sampler_res %>% 
    extract() %>%
    map(as.tibble) %>% 
    map2(names(.), 
         ., 
         ~set_names(.y, paste0(.x, '_', names(.y)))) %>%
    bind_cols
  
  c('muMutRNA',
    'muMutDNA',
    'muRefRNA',
    'muRefDNA') %>% 
    map(quant_norm_parameter, sample_df = sample_df) %>% 
    bind_cols %>% 
    transmute(transcriptional_shift = log(muMutRNA / muMutDNA) - log(muRefRNA / muRefDNA))
}

varInfo %<>% 
  mutate(transcriptional_shift_samples = mclapply(sampler_result, get_snp_TS_samples, mc.cores = 20),
         transcriptional_shift_HDI = map(transcriptional_shift_samples, ~pull(.x, transcriptional_shift) %>% mcmc %>% HPDinterval(prob = .99)),
         functional = map_lgl(transcriptional_shift_HDI, ~!between(0, .x[1], .x[2])),
         mean_transcriptional_shift = map_dbl(transcriptional_shift_samples, ~pull(.x, transcriptional_shift) %>% mean)) 

varInfo %>% 
  filter(functional) %>% 
  arrange(desc(abs(mean_transcriptional_shift))) %>% 
  filter(abs(mean_transcriptional_shift) < 10) %>% 
  save(file = '~/bayesianMPRA/outputs/functional_varInfo.RData')



# some gave ridiculous (10 to Inf) estimates. Others _changed sign_. I think these mostly caused by alleles/samples/barcodes where there are lots of 0's
