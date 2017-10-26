# Let's just try running a random forest with the histone marks on the Ulirsch
# data and looking at the variable importance metrics to see if we get ANTHING
# reasonable

library(caret)
library(tidyverse)
library(magrittr)

load("~/bayesianMPRA/analysis_data/varInfoWithHistoneMarkAnnotations.RData")

tmp = train(x = varInfo %>% dplyr::select(contains('K562')), y = varInfo %>% dplyr::select(transcription_shift) %>% pull)

varInfo %<>%
  mutate(fold = base::sample(1:10, 
                             size = nrow(.),
                             replace = TRUE))

train_histone_rf = function(fold_num){
  train_dat = varInfo %>% filter(fold != fold_num)
  test_dat = varInfo
}