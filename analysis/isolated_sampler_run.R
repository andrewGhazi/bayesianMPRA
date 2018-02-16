library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(fitdistrplus)
library(preprocessCore)
library(coda)
source('~/bayesianMPRA/analysis/mledistModified.R')

load('~/bayesianMPRA/analysis/outputs/objects_for_isolated_run.RData')

sampler_results = varInfo %>%
  group_by(construct) %>%
  nest %>%
  mutate(data = map2(data, construct, ~mutate(.x, construct = gsub('/', '-', gsub(' ', '_', .y))))) %$% 
  mclapply(data, run_sampler, mc.cores = 9)

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
    rstan::extract() %>%
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

load_construct_get_TS = function(construct){
  load(paste0('~/bayesianMPRA/outputs/ulirsch_sampler_results/', gsub('/', '-', gsub(' ', '_', construct)), '.RData'))
  
  get_snp_TS_samples(sampling_result)
}

varInfo %<>%
  mutate(transcriptional_shift_samples = map(construct, load_construct_get_TS),
         transcriptional_shift_HDI = map(transcriptional_shift_samples, ~pull(.x, transcriptional_shift) %>% mcmc %>% HPDinterval(prob = .99)),
         functional = map_lgl(transcriptional_shift_HDI, ~!between(0, .x[1], .x[2])),
         mean_transcriptional_shift = map_dbl(transcriptional_shift_samples, ~pull(.x, transcriptional_shift) %>% mean))

save(varInfo, file = '~/bayesianMPRA/outputs/varInfo_with_samples.RData')

bayesian_functional = varInfo %>%
  filter(functional) %>%
  arrange(desc(abs(mean_transcriptional_shift))) %>%
  filter(abs(mean_transcriptional_shift) < 10)

save(bayesian_functional, file = '~/bayesianMPRA/outputs/functional_varInfo.RData')



# some gave ridiculous (10 to Inf) estimates. Others _changed sign_. I think these mostly caused by alleles/samples/barcodes where there are lots of 0's
