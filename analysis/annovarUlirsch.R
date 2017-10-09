# Let's run annovar on the Ulirsch variants to try and find the best TS predictors

library(tidyverse)
library(magrittr)

load('outputs/varFuns.RData')

vcf = varFuns %>% select(chr, pos, ref, alt) %>% 
  unique %>% 
  rename(POS = pos) %>% 
  rename(CHROM = chr,
         REF = ref,
         ALT = alt) %>% 
  mutate(ID = '.',
         QUAL = '.', 
         FILTER = '.',
         INFO = '.') %>% 
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

vcf %>% 
  write_tsv('~/bayesianMPRA/outputs/ulirsch.vcf',
            col_names = FALSE) %>% 
  names %>% 
  paste0(collapse = '\t') %>% 
  paste0('#', ., '\\n') %>% 
  paste0('sed -i \'1s/^/', ., '/\' ~/bayesianMPRA/outputs/ulirsch.vcf')

ulirsch_anno = read_tsv('/usr/local/src/annovar-latest/ulirsch_dbnsfp_filter.hg19_dbnsfp30a_filtered',
                        col_names = FALSE,
                        na = '.') 

var

