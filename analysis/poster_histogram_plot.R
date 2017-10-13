library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(fitdistrplus)
library(preprocessCore)
library(ghzutils)

load('~/bayesianMPRA/outputs/objects_for_isolated_run.RData')

run_sampler = function(snp_data, n_iter = 3334){
  # snp_data - a data_frame with one row containing a column called countData and another called RNAgammaParams
  
  # Given a matrix of counts (rows = barcodes, columns = samples) and a
  # data_frame of by-allele-RNA-sample gamma hyperpriors, run the above Stan
  # model
  
  count_data = snp_data$countData[[1]]
  RNA_gamma_params = snp_data$RNAgammaParams[[1]]
  
  # Prepare count data matrices
  ref_DNA_mat = count_data %>% 
    filter(type == 'Ref') %>% 
    dplyr::select(-type) %>% 
    dplyr::select(contains('DNA')) %>% 
    as.matrix
  
  ref_RNA_mat = count_data %>% 
    filter(type == 'Ref') %>% 
    dplyr::select(-type) %>% 
    dplyr::select(contains('RNA')) %>% 
    as.matrix
  
  mut_DNA_mat  = count_data %>% 
    filter(type == 'Mut') %>% 
    dplyr::select(-type) %>% 
    dplyr::select(contains('DNA')) %>% 
    as.matrix
  
  mut_RNA_mat = count_data %>% 
    filter(type == 'Mut') %>% 
    dplyr::select(-type) %>% 
    dplyr::select(contains('RNA')) %>% 
    as.matrix
  
  # Prepare Gamma hyper-prior matrices
  mu_ref_RNA_hyper_params = RNA_gamma_params %>% 
    filter(type == 'Ref') %>% 
    pull(muGammaHyperPriors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  mu_mut_RNA_hyper_params = RNA_gamma_params %>% 
    filter(type == 'Mut') %>% 
    pull(muGammaHyperPriors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  phi_ref_RNA_hyper_params = RNA_gamma_params %>% 
    filter(type == 'Ref') %>% 
    pull(sizeGammaHyperPriors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  phi_mut_RNA_hyper_params = RNA_gamma_params %>% 
    filter(type == 'Mut') %>% 
    pull(sizeGammaHyperPriors) %>% 
    reduce(bind_rows) %>% 
    as.matrix() %>%
    t
  
  mu_DNA_hyper_params = margDNAPrior %>% 
    filter(grepl('mu', negBinParam)) %>% 
    dplyr::select(alphaEst, betaEst) %>% # OMG three pipes lining up
    as.matrix %>% 
    t
  
  phi_DNA_hyper_params = margDNAPrior %>% 
    filter(grepl('size', negBinParam)) %>% 
    dplyr::select(alphaEst, betaEst) %>% 
    as.matrix %>% 
    t
  
  # create input data list
  data_list = list(nRefBarcode = nrow(ref_DNA_mat),
                   nMutBarcode = nrow(mut_DNA_mat), 
                   nDNAblocks = ncol(ref_DNA_mat), 
                   nRNAblocks = ncol(ref_RNA_mat),
                   refDNAmat = ref_DNA_mat,
                   refRNAmat = ref_RNA_mat,
                   mutDNAmat = mut_DNA_mat,
                   mutRNAmat = mut_RNA_mat,
                   muRefRNAHyperParams = mu_ref_RNA_hyper_params,
                   phiRefRNAHyperParams = phi_ref_RNA_hyper_params,
                   muMutRNAHyperParams = mu_mut_RNA_hyper_params,
                   phiMutRNAHyperParams = phi_mut_RNA_hyper_params,
                   muDNAHyperParams = mu_DNA_hyper_params,
                   phiDNAHyperParams = phi_DNA_hyper_params)
  
  sampling(object = model,
           data = data_list,
           chains = 3,
           iter = n_iter,
           warmup = 500,
           thin = 1,
           verbose = FALSE) #friggin stan still verbose af
}

# get more samples
more_samples = varInfo %>% filter(construct == '9 136155000 2/3') %>% run_sampler(n_iter = 10000)

### Stan model -----
modelString_prior = "
data{
  int<lower=0> nRefBarcode ; // number of barcodes in ref allele
  int<lower=0> nMutBarcode ;
  int<lower=0> nDNAblocks ; // number of DNA replicates / blocks / samples / transfections
  int<lower=0> nRNAblocks ;
  int<lower=0> refDNAmat[nRefBarcode, nDNAblocks] ; // MPRA count matrix. Rows = barcodes, columns = samples/blocks
  int<lower=0> refRNAmat[nRefBarcode, nRNAblocks] ;
  int<lower=0> mutDNAmat[nMutBarcode, nDNAblocks] ;
  int<lower=0> mutRNAmat[nMutBarcode, nRNAblocks] ;
  real<lower=0> muRefRNAHyperParams[2, nRNAblocks] ; // gamma hyper-parameters on negative binomial parameters
  real<lower=0> phiRefRNAHyperParams[2, nRNAblocks] ;
  real<lower=0> muMutRNAHyperParams[2, nRNAblocks] ;
  real<lower=0> phiMutRNAHyperParams[2, nRNAblocks] ;
  real<lower=0> muDNAHyperParams[2, nDNAblocks] ;
  real<lower=0> phiDNAHyperParams[2, nDNAblocks] ;
}
parameters{
  real<lower=0> muRefDNA[nDNAblocks] ; //mean parameters for each block in each allele for each nucleic acid
  real<lower=0> muRefRNA[nRNAblocks] ;
  real<lower=0> muMutDNA[nDNAblocks] ;
  real<lower=0> muMutRNA[nRNAblocks] ;
  real<lower=0> phiRefDNA[nDNAblocks] ; // size parameters
  real<lower=0> phiRefRNA[nRNAblocks] ;
  real<lower=0> phiMutDNA[nDNAblocks] ;
  real<lower=0> phiMutRNA[nRNAblocks] ;
}
model{


 for (i in 1:nDNAblocks){

    // negative binomial parameters come from gamma hyper-priors
    muRefDNA[i] ~ gamma(muDNAHyperParams[1, i], muDNAHyperParams[2, i]) ; 
    phiRefDNA[i] ~ gamma(phiDNAHyperParams[1, i], phiDNAHyperParams[2, i]) ;
    muMutDNA[i] ~ gamma(muDNAHyperParams[1, i], muDNAHyperParams[2, i]) ; 
    phiMutDNA[i] ~ gamma(phiDNAHyperParams[1, i], phiDNAHyperParams[2, i]) ;

    // count data comes from the specified negative binomial
    // refDNAmat[,i] ~ neg_binomial_2(muRefDNA[i], phiRefDNA[i]) ;
    // mutDNAmat[,i] ~ neg_binomial_2(muMutDNA[i], phiMutDNA[i]) ;
  }

  for (i in 1:nRNAblocks){
    // negative binomial parameters come from gamma hyper-priors
    muRefRNA[i] ~ gamma(muRefRNAHyperParams[1, i], muRefRNAHyperParams[2, i]) ;
    muMutRNA[i] ~ gamma(muMutRNAHyperParams[1, i], muMutRNAHyperParams[2, i]) ;
    phiRefRNA[i] ~ gamma(phiRefRNAHyperParams[1, i], phiRefRNAHyperParams[2, i]) ;
    phiMutRNA[i] ~ gamma(phiMutRNAHyperParams[1, i], phiMutRNAHyperParams[2, i]) ;
    
    // count data comes from the specified negative binomial
    //refRNAmat[,i] ~ neg_binomial_2(muRefRNA[i], phiRefRNA[i]) ;
    //mutRNAmat[,i] ~ neg_binomial_2(muMutRNA[i], phiMutRNA[i]) ;
  }
}
"

prior_model = stan_model(model_code = modelString_prior)
prior_samples = sampling(object = prior_model,
                         data = data_list,
                         chains = 3,
                         iter = n_iter,
                         warmup = 500,
                         thin = 1,
                         verbose = FALSE)


more_samples = varInfo %>% filter(construct == '9 136155000 2/3') %>% run_sampler(n_iter = 10000)
#### Create the plot --------

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

p = get_snp_TS_samples(more_samples) %>% 
  ggplot(aes(transcriptional_shift)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 40, 
                 color = 'black', 
                 fill = 'grey60') +
  ggtitle('chr9:136155000 C to T from Ulirsch et al., Cell 2016') +
  xlab('Transcriptional Shift') +
  geom_segment(data = data_frame(x = .144761, xend = .872504, y = -.05, yend = -.05),
               inherit.aes = FALSE,
               mapping = aes(x = x, xend = xend, y = y, yend = yend),
               color = 'forestgreen', 
               lwd = 1.5) + 
  theme_pres()

ggsave(filename = '~/bayesianMPRA/outputs/plots/example_functional.png',
       plot = p,
       width = 10,
       height = 4.92)


