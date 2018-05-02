library(tidyverse)
library(magrittr)
library(rstan)

nb_reg_model = '
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=0> n_barcodes; // number in the WHOLE assay, not just per allele
  int<lower=0> rna_counts[n_barcodes, n_rna_samples]; 
  int<lower=0> dna_counts[n_barcodes, n_dna_samples]; 
  int<lower=0> n_alleles;
  int<lower=0, upper = n_alleles> allele[n_barcodes]; //allele indicator
  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> dna_depths[n_rna_samples];
} 
parameters {
  real<lower=0> rna_m_a;
  real<lower=0> rna_m_b;
  real<lower=0> dna_m_a;
  real<lower=0> dna_m_b;
  real<lower=0> rna_p_a;
  real<lower=0> rna_p_b;
  real<lower=0> dna_p_a;
  real<lower=0> dna_p_b;

  real r_m_i[n_alleles,n_rna_samples];
  real r_p_i[n_alleles,n_rna_samples];
  real d_m_i[n_alleles,n_rna_samples];
  real d_p_i[n_alleles,n_rna_samples];
}
model {

  
  rna_m_a ~ gamma(20,20);
  rna_m_b ~ gamma(20,20);
  rna_p_a ~ gamma(20,20);
  rna_p_b ~ gamma(20,20);

  dna_m_a ~ gamma(20,20);
  dna_m_b ~ gamma(20,20);
  dna_p_a ~ gamma(20,20);
  dna_p_b ~ gamma(20,20);

  for (s in 1:n_rna_samples) {
    r_m_i[allele,s] ~ gamma(rna_m_a, rna_m_b);
    r_p_i[allele,s] ~ gamma(rna_p_a, rna_p_b);
    rna_counts[allele, s] ~ neg_binomial_2(r_m_i[allele], r_p_i[allele]);
  }

  for (s in 1:n_dna_samples){
    d_m_i[allele,s] ~ gamma(dna_m_a, dna_m_b);
    d_p_i[allele,s] ~ gamma(dna_p_a, dna_p_b);
    dna_counts[allele, s] ~ neg_binomial_2(d_m_i[allele], d_p_i[allele]);
  }
}
'

nb_reg_object = stan_model(model_code = nb_reg_model)

data_list = list()

## basal mean ----

nb_basal_model = '
data {
  int<lower=0> n_rna_samples;
  int<lower=0> n_dna_samples;
  int<lower=0> n_barcodes; // number in the WHOLE assay, not just per allele
  int<lower=0> rna_counts[n_barcodes, n_rna_samples]; 
  int<lower=0> dna_counts[n_barcodes, n_dna_samples]; 
  int<lower=0> n_alleles;
  int<lower=0, upper = n_alleles> allele[n_barcodes]; //allele indicator
  real<lower=0> rna_depths[n_rna_samples]; 
  real<lower=0> dna_depths[n_rna_samples];
} 
parameters {
  real<lower=0> rna_m_a;
  real<lower=0> rna_m_b;
  real<lower=0> dna_m_a;
  real<lower=0> dna_m_b;
  real<lower=0> rna_p_a;
  real<lower=0> rna_p_b;
  real<lower=0> dna_p_a;
  real<lower=0> dna_p_b;

  real r_m_i[n_alleles,n_rna_samples];
  real r_p_i[n_alleles,n_rna_samples];
  real d_m_i[n_alleles,n_rna_samples];
  real d_p_i[n_alleles,n_rna_samples];
}
model {


  rna_m_a ~ gamma(20,20);
  rna_m_b ~ gamma(20,20);
  rna_p_a ~ gamma(20,20);
  rna_p_b ~ gamma(20,20);
  
  dna_m_a ~ gamma(20,20);
  dna_m_b ~ gamma(20,20);
  dna_p_a ~ gamma(20,20);
  dna_p_b ~ gamma(20,20);
  
  for (s in 1:n_rna_samples) {
    r_m_i[allele,s] ~ gamma(rna_m_a, rna_m_b);
    r_p_i[allele,s] ~ gamma(rna_p_a, rna_p_b);
    rna_counts[allele, s] ~ neg_binomial_2(r_m_i[allele], r_p_i[allele]);
  }

  for (s in 1:n_dna_samples){
    d_m_i[allele,s] ~ gamma(dna_m_a, dna_m_b);
    d_p_i[allele,s] ~ gamma(dna_p_a, dna_p_b);
    dna_counts[allele, s] ~ neg_binomial_2(d_m_i[allele], d_p_i[allele]);
  }
}
'