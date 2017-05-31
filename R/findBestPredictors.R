varFuns %>% 
  filter(abs(TS) > .33, abs(TS) < 2) %>%
  select(TS, eigen:DeepSeaDnaase, BroadK562H3k4me1:SydhK562H3k27me3b) %>%
  gather(predictor, val, -TS) %>% 
  na.omit %>% 
  group_by(predictor) %>% 
  summarise(corr = cor(TS, val)) %>% 
  arrange(desc(abs(corr)))

getRefNNdiffs = function(NNs, meanRefMu) {
  meanRefMu - varFuns[NNs,]$meanRefMu
}

getRefRandDiffs = function(meanRefMu) {
  meanRefMu - varFuns[sample(1:7803, 30),]$meanRefMu
}

getMutNNdiffs = function(NNs, meanMutMu) {
  meanMutMu - varFuns[NNs,]$meanMutMu
}

getMutRandDiffs = function(meanMutMu) {
  meanMutMu - varFuns[sample(1:7803, 30),]$meanMutMu
}

varFuns %<>% 
  mutate(RefNNdiffs = map2(kNN, meanRefMu, getRefNNdiffs),
         RefRandDiffs = map(meanRefMu, getRefRandDiffs),
         MutNNdiffs = map2(kNN, meanMutMu, getMutNNdiffs),
         MutRandDiffs = map(meanMutMu, getMutRandDiffs))

tmp = varFuns %>% 
  filter(abs(transcriptionalShift) > .33) %>% 
  select(contains('iffs')) %>% 
  unnest %>% 
  gather(src, diffs)

tmp %>% 
  ggplot(aes(diffs)) + 
  geom_histogram(bins = 50) + 
  facet_grid(src ~ .)

tmp %>% 
  group_by(src) %>% 
  summarise(sd = sd(diffs),
            mn = mean(diffs))
