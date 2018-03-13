## Let's annotate the Tewhey data

library(parallel)
library(tidyverse)


dir = '/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/'

file_names = list.files(dir, 
                        pattern = '.expanded$') %>% 
  grep(pattern = 'ctrl|HepG2', x = ., value = TRUE)

tewhey_dat = mclapply(file_names,
                         function(t_file){
                           read_tsv(paste0(dir, t_file), 
                                    col_names = c('snp_allele', 'barcode', 'count')) %>% 
                             mutate(sample = t_file) %>% 
                             separate(snp_allele, into = c('snp_id', 'allele'), sep = stringr::regex('(?=[AB]$)')) %>% 
                             mutate(snp_id = gsub('_$', '', snp_id))
                         },
                         mc.cores = 11) %>% 
  bind_rows

save(tewhey_subset, file = '~/bayesianMPRA/analysis_data/tewhey_subset.RData')


