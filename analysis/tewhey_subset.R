
dir = '/mnt/bigData2/andrew/MPRA/Tewhey/indivTags/'

file_names = list.files(dir, 
                        pattern = '.expanded$') %>% 
  grep(pattern = 'ctrl|HepG2', x = ., value = TRUE)

tewhey_subset = mclapply(file_names,
                 function(t_file){
                   read_tsv(paste0(dir, t_file), 
                            col_names = c('snp_allele', 'barcode', 'count')) %>% 
                     mutate(sample = t_file) %>% 
                     separate(snp_allele, into = c('snp_id', 'allele'), sep = stringr::regex('(?=[AB]$)')) %>% 
                     mutate(snp_id = gsub('_$', '', snp_id))
                 },
                 mc.cores = 11) %>% 
  bind_rows %>% 
  filter(snp_id %in% base::sample(unique(snp_id), size = 2000)) %>% #randomly choose 5000 alleles as a subset
  unique %>% # there are infrequently some duplicate rows in the tewhey data
  mutate(sample = gsub('Geuv_90K_', '', sample) %>% gsub('.tag.ct.indiv.expanded', '', .))

save(tewhey_subset, file = '~/bayesianMPRA/analysis_data/tewhey_subset.RData')

