# Let's simulate data based on simulated distributions of effect sizes and see
# how the power shakes out

# TODO resample nb mle estimates, rather than TS. That would probably be a more realistic simulation

library(tidyverse)

n_allele = 1000
n_dna = 3
n_rna = 5
n_barcode = 20

load("~/bayesianMPRA/analysis_data/varInfoWithHistoneMarkAnnotations.RData")
ts_dens = density(varInfo$transcription_shift, n = 1024)

sample_depths = data_frame(sample = c(paste0('dna', 1:n_dna), paste0('rna', 1:n_rna)),
                           depth = floor(runif(n_rna+n_dna, 4e5, 2e6)))

sample_ts_dens = function(n){
  sample(ts_dens$x,
         size = n,
         replace = TRUE,
         prob = ts_dens$y)
}

simulate_count = function(true_shift, n_dna, n_rna, n_barcode){
  data_frame(allele = rep(c('ref', 'mut'), each = n_barcode*(n_dna+n_rna)),
             sample = rep(c(paste0('dna', 1:n_dna), paste0('rna', 1:n_rna)), times = 2*n_barcode))
}

simulate_counts = function(shift_vec, n_dna, n_rna, n_barcode){
  
}

sim_data = data_frame(allele = map_chr(1:n_allele,
                                       ~paste0(base::sample(letters, size = 10, replace = TRUE), collapse = '')),
                      true_ts = sample_ts_dens(n_allele),
                      sim_counts = simulate_counts(true_ts))
