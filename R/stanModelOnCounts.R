library(tidyverse)
library(magrittr)
library(rstan)


modelString = "
data{
  int<lower=0> nRefBarcode ;
  int<lower=0> nMutBarcode ;
  int<lower=0> nDNAblocks ;
  int<lower=0> nRNAblocks ;
  int<lower=0> refDNAmat[nRefBarcode, nDNAblocks] ;
  int<lower=0> refRNAmat[nRefBarcode, nRNAblocks] ;
  int<lower=0> mutDNAmat[nMutBarcode, nDNAblocks] ;
  int<lower=0> mutRNAmat[nMutBarcode, nRNAblocks] ;
  real<lower=0> muHyperParams[2] ;
  real<lower=0> phiHyperParams[2] ;
}
parameters{
  real<lower=0> muRefDNA[nDNAblocks] ;
  real<lower=0> muRefRNA[nRNAblocks] ;
  real<lower=0> muMutDNA[nDNAblocks] ;
  real<lower=0> muMutRNA[nRNAblocks] ;
  real<lower=0> phiRefDNA[nDNAblocks] ;
  real<lower=0> phiRefRNA[nRNAblocks] ;
  real<lower=0> phiMutDNA[nDNAblocks] ;
  real<lower=0> phiMutRNA[nRNAblocks] ;
}
model{
  
  muRefDNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;
  muRefRNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;
  muMutDNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;
  muMutRNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;

  phiRefDNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
  phiRefRNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
  phiMutDNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
  phiMutRNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
  
  for (i in 1:nDNAblocks){
    refDNAmat[,i] ~ neg_binomial_2(muRefDNA[i], phiRefDNA[i]) ;
    mutDNAmat[,i] ~ neg_binomial_2(muMutDNA[i], phiMutDNA[i]) ;
  }
  for (i in 1:nRNAblocks){
    refRNAmat[,i] ~ neg_binomial_2(muRefRNA[i], phiRefRNA[i]) ;
    mutRNAmat[,i] ~ neg_binomial_2(muMutRNA[i], phiMutRNA[i]) ;
  }
}
"

model = stan_model(model_code = modelString)

set.seed(1280)

dir = "/mnt/labhome/andrew/MPRA/paper_data/"

exampleData = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                         delim = "\t",
                         col_names = T,
                         col_types = cols(chr = "c")) %>% 
  dplyr::select(matches('construct|CTRL|DNA|type')) %>% 
  filter(construct == '1 155271258 1/3') %>% 
  mutate(bcIndex = 1:nrow(.)) 

dataList = list(nRefBarcode = sum(exampleData$type == 'Ref'),
                nMutBarcode = sum(exampleData$type == 'Mut'),
                nDNAblocks = colnames(exampleData) %>% grepl('DNA', .) %>% sum,
                nRNAblocks = colnames(exampleData) %>% grepl('RNA', .) %>% sum,
                refDNAmat = exampleData %>% filter(type == 'Ref') %>% dplyr::select(contains('DNA')) %>% as.data.frame %>% as.matrix,
                refRNAmat = exampleData %>% filter(type == 'Ref') %>% dplyr::select(contains('RNA')) %>% as.data.frame %>% as.matrix,
                mutDNAmat = exampleData %>% filter(type == 'Mut') %>% dplyr::select(contains('DNA')) %>% as.data.frame %>% as.matrix,
                mutRNAmat = exampleData %>% filter(type == 'Mut') %>% dplyr::select(contains('RNA')) %>% as.data.frame %>% as.matrix,
                muHyperParams = c(.01, .01),
                phiHyperParams = c(.01, .01))

stanFit = sampling(object = model,
                   data = dataList,
                   chains = 3, 
                   iter = 3334,
                   warmup = 500,
                   thin = 1)

sampleList = rstan::extract(stanFit)
paramNames = names(sampleList) %>% .[-length(.)]

paramSamplesToDF = function(paramName){
  sampleList[[paramName]] %>% 
    as.data.frame %>% 
    purrr::set_names(., nm = paste0(paramName, 'Block', 1:ncol(.))) %>% 
    as.tbl
}

allSamples = paramNames %>% 
  map(paramSamplesToDF) %>% 
  purrr::reduce(bind_cols)

meanParameterSamples = paramNames %>% 
  map(~allSamples %>% select(contains(.x)) %>% rowMeans %>% as.data.frame %>% set_names(.x) %>% as.tbl) %>% #take the mean within blocks
  purrr::reduce(bind_cols)

meanParameterSamples %<>%
  mutate(TS = log(muMutRNA/muMutDNA) - log(muRefRNA/muRefDNA))
