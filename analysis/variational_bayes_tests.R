run_vb = function(snp_data, n_iter = 3334){
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
  
  sampling_result = sampling(object = model,
                             data = data_list,
                             chains = 3,
                             iter = n_iter,
                             warmup = 500,
                             thin = 1,
                             verbose = FALSE) #friggin stan still verbose af
  
  save(list = c('sampling_result' , 'snp_data'),
       file = paste0('~/bayesianMPRA/outputs/ulirsch_sampler_results/', snp_data$construct, '.RData'))
  
  return('finished')
}