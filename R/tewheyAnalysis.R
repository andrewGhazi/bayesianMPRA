library(tidyverse)
library(readr)
library(stringr)

tewhey = read_tsv('/mnt/bigData2/andrew/MPRA/Tewhey/GSE75661_79k_collapsed_counts.txt') %>% 
  separate(col = Oligo, into = c('oligo', 'allele'), sep = '(?=[AB]$)')

# Check number of 0's in each column
tewhey %>% map_int(~sum(.x == 0))

DNAmeans = tewhey %>% 
  gather(transfection, count, -oligo, -allele) %>% 
  group_by(transfection) %>% 
  mutate(count = 1e6*count / sum(count)) %>% #Normalize by transfection depth
  ungroup %>% group_by(oligo, allele) %>% summarise(meanPlasmid = mean(count[grepl('Plasmid', transfection)])) %>% ungroup
