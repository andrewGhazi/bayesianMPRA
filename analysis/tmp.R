tmp = metrop(varFuns$postFun[[981]], initial = list(refMu = 0, mutMu = 0, sig = 1), nbatch = 3, blen = 1000, nspac = 5)
# tmp = expand.grid(seq(-10,10, by = .1), seq(-10,10, by = .1), seq(.1, 5, by = .1)) %>% as.tbl
# library(parallel)
# tmp$postDens = mclapply(seq_along(tmp$Var1),
#                         function(x) {
#                           varFuns$postFun[[981]](tmp$Var1[x], tmp$Var2[x], tmp$Var3[x])
#                         }, mc.cores = 12) %>% unlist

tmpfun = function(x) {
  return(x**2 + (x - 2)**2)
}

res$batch %>% 
  as.data.frame %>% 
  as.tbl %>% 
  set_names(c('refMu', 'mutMu', 'sig')) %>% 
  mutate(TS = mutMu - refMu) %>% 
  ggplot(aes(x = TS, y = ..density..)) + 
  geom_histogram(bins = 50)
