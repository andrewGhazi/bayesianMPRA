#Let's compute the pooled standard deviation
load('data/varInfoWithHistoneMarkAnnotations.RData')

library(purrr)

getConstructRefStdDev = function(constr){
  MPRA.qnactivity %>% filter(construct == constr, type == 'Ref') %>% .$qnact %>% sd
}

varInfo %<>% 
  mutate(refActStdDev = map_dbl(construct, getConstructRefStdDev))

save(varInfo, file = 'data/varInfoWithHistoneMarkAnnotations.RData')