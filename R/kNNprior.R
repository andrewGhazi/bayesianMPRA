library(tidyverse)
library(magrittr)
library(stringr)

load('data/varInfoWithHistoneMarkAnnotations.RData')
names(varInfo)[12] = 'transcriptionalShift'

preds = varInfo %>% select(contains('Broad'), contains('Sydh')) %>% 
  gather(pred, val) %>% 
  group_by(pred) %>% 
  mutate(val = scale(map_dbl(val, ~ifelse(.x == 0, .1*min(val[val>0]), .x)))) %>%  #replace 0's with one tenth the minimal value for the log
  ungroup
  

