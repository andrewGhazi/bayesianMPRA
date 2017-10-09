library(tidyverse)
library(magrittr)
source('~/utility/scripts/utilityFunctions.R')

load('data/varInfoWithHistoneMarkAnnotations.RData')

read_regulome_dat = function(chr_num){ # this doesn't work
  reg_dat = read_tsv(paste0('data/RegulomeDB/chr', chr_num),
                     col_names = FALSE)
  
  reg_dat %<>%
    select(1:9, ncol(reg_dat))
}