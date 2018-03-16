## Let's annotate the Tewhey data

library(parallel)
library(data.table)
library(magrittr)
library(stringr)
library(tidyverse)

tewhey_snps = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/GSE75661_79k_collapsed_counts.txt', # this data is already summed across barcodes, I'm just using it to get the sample depths
                       col_names = TRUE) %>% 
  select(Oligo) %>% 
  separate(Oligo, into = c('snp_id', 'allele'), sep = stringr::regex('(?=[AB]$)')) %>%
  mutate(snp_id = gsub(pattern = '_$', replacement = '', x = snp_id, perl = TRUE)) %>% 
  filter(grepl('rs', snp_id)) %>%  # Most of those removed are indels (eg chr1:150824527:I_RC, 5452 of them) or formats I don't understand (chr6:32629802,  MERGED_DEL_2_2432_RC, 102 of them)
  select(-allele) %>% 
  mutate(rs_id = str_extract(snp_id, 'rs[0-9]+$')) %>% 
  select(rs_id) %>% 
  unique %>% 
  na.omit

write_tsv(tewhey_snps,
          path = 'analysis_outputs/tewhey_snps.tsv',
          col_names = FALSE)

# Command to get the appropriate annovar input of the Tewhey rsid's. convert2annovar.pl doesn't work, so the github link below was adapted to this command.
# grep -w -f ~/bayesianMPRA/analysis_outputs/tewhey_snps.tsv hg19_avsnp147.txt > ~/bayesianMPRA/analysis_outputs/tewhey_snps.avinput
# https://github.com/WGLab/doc-ANNOVAR/issues/16

# Annotate these with annovar

tewhey_snps %>%
  select(snp_id) %>%
  unique %>% 
  mutate(snp_id = gsub('_RC$', '', snp_id, perl = TRUE)) %>% # I believe these refer to reverse complements 
  unique

#### Depth normalization ----
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
