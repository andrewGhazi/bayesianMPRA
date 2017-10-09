# Need fasta's to get quotes for plasmids for luciferase validation

library(tidyverse)
library(stringr)
library(magrittr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqinr)

genome = BSgenome.Hsapiens.UCSC.hg19

bfunc = read_tsv('~/bayesianMPRA/outputs/bayesian_functional_variants.tsv')

to_test = bfunc %>% 
  arrange(priority) %>% 
  .[1:3,] %>% 
  mutate(reverse_gene = c(FALSE, TRUE, FALSE),
         variant_pos = str_extract(construct, '[0-9]/[0-9]'))

context_from_dat = function(construct_dat){
  
  if (construct_dat$variant_pos == '1/3') {
    
    start_pos = construct_dat$pos - 48
    end_pos = construct_dat$pos + 96
    
    if (construct_dat$reverse_gene) {
      start_pos = construct_dat$pos - 96
      end_pos = construct_dat$pos + 48
    }
    
  } else if (construct_dat$variant_pos == '2/3') {
    
    start_pos = construct_dat$pos - 96
    end_pos = construct_dat$pos + 48
    
    if (construct_dat$reverse_gene) {
      start_pos = construct_dat$pos - 48
      end_pos = construct_dat$pos + 96
    }
    
  } else if (construct_dat$variant_pos == '1/2') { 
    
    start_pos = construct_dat$pos - 72
    end_pos = construct_dat$pos + 72
    
  } else {
    
    stop('incorrect position given')
    
  }
  
  construct_dat %<>% 
    mutate(ref_seq = subseq(genome[[paste0('chr', construct_dat$chr)]],
                                  start = start_pos,
                                  end = end_pos) %>% toString)
  
  pos_to_replace = 144 * map_dbl(strsplit(construct_dat$variant_pos, '/'), ~as.numeric(.x[[1]]) / as.numeric(.x[[2]])) + 1
  
  if (construct_dat$reverse_gene) {
    pos_to_replace = 144 * (1 - map_dbl(strsplit(construct_dat$variant_pos, '/'), ~as.numeric(.x[[1]]) / as.numeric(.x[[2]]))) + 1
  }
  
  construct_dat %<>%
    mutate(mut_seq = str_c(str_sub(construct_dat$ref_seq, 1, pos_to_replace - 1),
                           construct_dat$alt,
                           str_sub(construct_dat$ref_seq, pos_to_replace + 1, nchar(construct_dat$ref_seq))))
  
  if (construct_dat$reverse_gene) {
    construct_dat %<>%
      mutate(ref_seq = ref_seq %>% DNAString() %>% reverseComplement() %>% toString,
             mut_seq = mut_seq %>% DNAString() %>% reverseComplement() %>% toString)
  }
  
  return(construct_dat)
}

seqs_to_write = to_test %>% 
  group_by(construct) %>% 
  nest %>% 
  mutate(context_dat = map(data, context_from_dat)) %>% 
  dplyr::select(-data) %>% 
  unnest %>% 
  select(construct, ref_seq, mut_seq) %>% 
  gather(seq_type, seq, -construct) %>% 
  mutate(construct = gsub(' ', '_', construct), 
         seq = map()) %>% 
  unite(seq_name, c('construct', 'seq_type'))

write.fasta(sequences = seqs_to_write$seq %>% str_extract_all('[ACGT]') %>% map(tolower),
            names = seqs_to_write$seq_name,
            file.out = '~/bayesianMPRA/outputs/luc_sequences.fasta')

