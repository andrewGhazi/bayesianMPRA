library(tidyverse)
library(stringr)
library(magrittr)

load('~/bayesianMPRA/outputs/varFuns.RData')
load('~/bayesianMPRA/data/gatheredMPRA.RData')

getWilcoxP = function(constr) {
  MPRA.qnactivity %>% 
    filter(construct == constr) %>% 
    summarise(pVal = wilcox.test(qnact[type == 'Ref'], qnact[type == 'Mut'])$p.value) %>% 
    .$pVal
}

varFuns %<>% 
  mutate(wilcoxP = map_dbl(construct, getWilcoxP))


save(varFuns, file = '~/bayesianMPRA/outputs/varFuns.RData')
