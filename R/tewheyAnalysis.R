library(tidyverse)
library(readr)
library(stringr)


tewhey = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/GSE75661_79k_collapsed_counts.txt') %>% 
  separate(col = Oligo, into = c('oligo', 'allele'), sep = '(?=[AB]$)') %>% 
  mutate(rs = str_extract(oligo, 'rs[0-9]+'))

# Check number of 0's in each column
tewhey %>% map_int(~sum(.x == 0))

tewhey %>% select(rs) %>% unique %>% write_tsv('~/bayesianMPRA/outputs/tewheyRS.txt', col_names = FALSE)

DNAmeans = tewhey %>% 
  gather(transfection, count, -oligo, -allele) %>% 
  group_by(transfection) %>% 
  mutate(count = 1e6*count / sum(count)) %>% #Normalize by transfection depth
  ungroup %>% group_by(oligo, allele) %>% summarise(meanPlasmid = mean(count[grepl('Plasmid', transfection)])) %>% ungroup

#Recapitulating counts
tewheyExample = tmp = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/Geuv_90K_NA12878.r1.tag.ct.indiv.expanded', col_names = c('allele', 'barcode', 'count'))

tewheyExample %>% 
  head %>% 
  mutate(myCount = map_int(barcode, ~as.integer(system(paste0('grep -c ', 
                                                              .x %>% DNAString %>% reverseComplement() %>% toString, 
                                                              ' /mnt/bigData2/andrew/MPRA/Tewhey/SRR2972757.fastq'), 
                                                       intern = TRUE))))
