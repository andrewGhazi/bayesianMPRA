library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(fitdistrplus)

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

more_samples = varInfo %>% filter(construct == '9 136155000 2/3') %>% run_sampler(n_iter = 10000)

p = get_snp_TS_samples(tmp) %>% 
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
