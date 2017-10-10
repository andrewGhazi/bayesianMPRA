
library(tidyverse)
library(biomaRt)
library(stringr)
select = dplyr::select

snpmart = useMart("ENSEMBL_MART_SNP",  
                  host="grch37.ensembl.org",
                  dataset = 'hsapiens_snp')

load('~/designMPRA/outputs/tewheyDepthAdjDNAcounts.RData')
load('~/designMPRA/outputs/tewheyDNAcounts.RData')

#### Prepare to annotate with DeepSea --------
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

deepsea_output = read_csv('~/bayesianMPRA/data/tewhey_deepsea/infile.vcf.out.logfoldchange') %>% 
  select(2:6, contains('HepG2')) %>% # Let's just try the HepG2 samples
  select(2:6, contains('DNase'))

#### Prepare for BMPRA -------
read_rna_fun = function(rna_file){
  read_tsv(paste0('/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/', rna_file),
           col_names = c('allele', 'barcode', 'count')) %>% 
    mutate(file = gsub('Geuv_90K_|.tag.ct.indiv.expanded', '', rna_file))
}

rna_counts = list.files('/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/',
                        pattern = 'HepG2') %>% 
  map_dfr(read_rna_fun)

full_dat = rna_counts %>% 
  bind_rows(dnaCounts %>% rename(file = src)) %>% 
  mutate(type = str_extract(allele, 'A$|B$'),
         snp = gsub('A$|B$','', allele) %>% gsub('_$', '', .)) %>% 
  select(-allele) %>% 
  group_by(snp) %>% 
  nest %>% 
  mutate(rs_id = str_extract(snp, 'rs[0-9]+')) %>% 
  left_join(deepsea_output, by = c('rs_id' = 'name'))

spread_and_omit = function(tidy_counts){
  # Need to put the counts into a matrix and throw out those that aren't present
  # in all samples
  
  tidy_counts %>% 
    unique %>% # Tewhey's data has some duplicate entries sometimes
    group_by(type) %>% 
    spread(file, count) %>% 
    ungroup %>% 
    na.omit 
    
}


