## Let's annotate the Tewhey data

library(parallel)
library(data.table)
library(magrittr)
library(tidyverse)

sample_sums = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/GSE75661_79k_collapsed_counts.txt', # this data is already summed across barcodes, I'm just using it to get the sample depths
                      col_names = TRUE) %>%
  separate(Oligo, into = c('snp_id', 'allele'), sep = stringr::regex('(?=[AB]$)')) %>%
  mutate(snp_id = gsub(pattern = '_$', replacement = '', x = snp_id, perl = TRUE)) %>%
  gather(sample, count, -snp_id, -allele)

tewhey_depths = sample_sums %>%
  group_by(sample) %>%
  summarise(depth = sum(count))

dna_depths = tewhey_depths %>% 
  filter(grepl('Plasmid', sample)) %>% 
  mutate(sample = str_extract(sample, 'r[0-9]')) %>% 
  rename(file = sample) %>% 
  as.data.table %>% 
  .[,file := factor(file)] %>% 
  setkey(file)

dir = '/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/'
dna_data_files = list.files(dir,
                            pattern = 'ctrl')

dna_counts = map(dna_data_files, ~fread(paste0(dir, .x),
                                        col.names = c('snp_id', 'barcode', 'count'))[,file := str_extract(.x, 'r[0-9]')]) %>% 
  bind_rows %>% 
  .[,file := factor(file)] %>% 
  setkey(file)

dna_norm_depths = dna_counts[dna_depths] %>% 
  .[,norm_depth := count * 1e6 / depth] %>% 
  .[,.(snp_id = snp_id[1], mean_norm_depth = mean(norm_depth)), by = barcode]

dna_norm_depths %>% ggplot(aes(mean_norm_depth)) + geom_histogram(bins = 50) + scale_x_log10() + geom_vline(xintercept = .03, lty = 2)

#### This data is already summed across barcodes
# tewhey_dat = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/GSE75661_79k_collapsed_counts.txt',
#                       col_names = TRUE) %>%
#   separate(Oligo, into = c('snp_id', 'allele'), sep = stringr::regex('(?=[AB]$)')) %>%
#   mutate(snp_id = gsub(pattern = '_$', replacement = '', x = snp_id, perl = TRUE)) %>%
#   gather(sample, count, -snp_id, -allele)
# 
# tewhey_depths = tewhey_dat %>%
#   group_by(sample) %>%
#   summarise(depth = sum(count))
# 
# tewhey_dna = tewhey_dat %>% 
#   filter(grepl('Plasmid', sample))
# 
# tehey_dna %>% 
#   filter(grepl('Plasmid', sample)) %>% 
#   ggplot(aes(count)) +
#   geom_histogram(bins = 40) + 
#   facet_wrap('sample') +
#   scale_x_log10()
# 
