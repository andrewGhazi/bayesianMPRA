
library(tidyverse)
library(biomaRt)
library(stringr)

snpmart = useMart("ENSEMBL_MART_SNP",  
                  host="grch37.ensembl.org",
                  dataset = 'hsapiens_snp')

load('~/designMPRA/outputs/tewheyDNAcounts.RData')

select = dplyr::select
tewhey_alleles = dnaCounts %>% 
  select(allele) %>% 
  unique %>% 
  mutate(rs_id = str_extract(allele, 'rs[0-9]+'))

# chr_pos_given = tewhey_alleles %>% # No easy way to figure out the alleles from these
#   filter(is.na(rs_id)) %>% 
#   select(-rs_id) %>% 
#   mutate(chr = str_extract(allele, 'chr[0-9]+') %>% str_replace('chr', '') %>% as.integer,
#          pos = str_extract(allele, '[0-9]{3,}') %>% as.integer)

rs_alleles = tewhey_alleles %>% 
  filter(!is.na(rs_id)) %>% 
  select(rs_id) %>% 
  unique
  
rs_to_pos = getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
      filters = 'snp_filter',
      values = rs_alleles %>% pull(rs_id),
      mart = snpmart) %>% 
  as.tibble

rs_to_annotate = rs_to_pos %>% 
  filter(chrom_start == chrom_end) %>% # this only cuts out 13
  rename(CHROM = chr_name,
         POS = chrom_start) %>% 
  separate(allele, into = c('REF', 'ALT'), sep = '/', extra = 'drop') # dropping third alleles for now

rs_to_vcf = rs_to_annotate %>% 
  rename(ID = refsnp_id) %>% 
  select(CHROM, POS, ID, REF, ALT) %>% 
  filter(ALT %in% c('A', 'C', 'G', 'T')) %T>%
  write_tsv(path = '~/bayesianMPRA/outputs/tewhey.vcf',
            col_names = TRUE) # Add the # in the actual file manually


