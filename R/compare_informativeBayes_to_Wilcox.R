library(ggbeeswarm)
select = dplyr::select

# wilcox
dir = "/mnt/labhome/andrew/MPRA/paper_data/"

UMPRA = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                   delim = "\t",
                   col_names = T,
                   col_types = cols(chr = "c"))

depthNormalize = function(sampleCounts){
  sampleCounts*1e6/sum(sampleCounts)
}

ulirsch_activities = UMPRA %>% 
  mutate_at(vars(contains('K562')), depthNormalize) %>% 
  mutate(dnaMean = (K562_minP_DNA1 + K562_minP_DNA2)/2) %>% 
  filter(dnaMean > .13) %>% 
  select(-K562_minP_DNA1, -K562_minP_DNA2) %>% 
  select(chr:K562_CTRL_minP_RNA6, dnaMean) %>% 
  gather(sample, depthAdjCount, K562_CTRL_minP_RNA1:K562_CTRL_minP_RNA6) %>% 
  mutate(activity = log(depthAdjCount / dnaMean))

ulirsch_control_wilcox = ulirsch_activities %>% 
  group_by(construct) %>% 
  nest %>% 
  mutate(has_both_alleles = map_lgl(data, ~(sum(.x$type == 'Ref') > 1 & sum(.x$type == 'Mut') > 1))) %>% 
  filter(has_both_alleles) %>% 
  mutate(wilcox_test = map(data, ~wilcox.test(.x$activity[.x$type == 'Ref'], .x$activity[.x$type == 'Mut'], verbose = FALSE)),
         wilcox_p = map_dbl(wilcox_test, ~.x$p.value),
         wilcox_q = p.adjust(wilcox_p, method = 'fdr'),
         wilcox_functional = wilcox_q < .01) %>% 
  arrange(wilcox_q)
  
ulirsch_control_wilcox %>% filter(wilcox_functional)

detected_only_by_bayes = varInfo %>% 
  filter(functional) %>% 
  arrange(desc(abs(mean_transcriptional_shift))) %>% 
  filter(abs(mean_transcriptional_shift) < 10) %>% 
  left_join(ulirsch_control_wilcox, by = 'construct') %>%
  select(-(gkmerDelta:DeepSeaDnaase), 
         -(weights:transcriptional_shift_samples),
         -(has_both_alleles:wilcox_test)) %>% 
  filter(functional, !wilcox_functional) %>% 
  arrange(abs(mean_transcriptional_shift))

detected_only_by_bayes %>% 
ggplot(aes(mean_transcriptional_shift)) + geom_rug() + geom_histogram()

ulirsch_activities %>% 
  filter(construct %in% detected_only_by_bayes$construct) %>% 
  ggplot(aes(type, activity)) + 
  geom_jitter(height = 0, size = .75, alpha = .3, width = .25) +
  facet_wrap('construct')

priority_vec = c(1,1)
priority = 1
while(priority_vec < 47){
  priority = priority + 1
  new_priority_vec = vector(mode = 'integer', length = 2*length(priority_vec) - 1)
  
  for(i in 1:priority)
}

detected_only_by_bayes %>% 
  mutate(priority = c(1,7,6,7,5,7,6,7,4,7,6,7,5,7,6,7,3,7,6,7,5,7,6,7,4,7,6,7,5,6,2,6,5,6,4,6,5,6,3,6,5,6,4,6,5,6,1)) %>% # lol
  select(-(data)) %>% 
  arrange(priority) %>% 
  mutate(HDI_lower = map_dbl(transcriptional_shift_HDI, ~.x[1]),
         HDI_upper = map_dbl(transcriptional_shift_HDI, ~.x[2])) %>% 
  select(-transcriptional_shift_HDI) %>% 
  rename(traditional_transcription_shift = transcriptionalShift) %>% 
  write_tsv('~/bayesianMPRA/outputs/bayesian_functional_variants.tsv')

