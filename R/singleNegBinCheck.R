# Let's do a back of the envelope check on how well the nbinom works for counts of a single allele in a single sample

# > pltMPRA %>% filter(rs == 'ALAS21') %>% spread(sample, count) %>% dplyr::select(type, contains('run')) %>% filter(type == 'Ref', run1_EB9223_GCCAAT > 10) %>% pull(run1_EB9215_GATCAG) %>% mean
# [1] 86.59375
# > pltMPRA %>% filter(rs == 'ALAS21') %>% spread(sample, count) %>% dplyr::select(type, contains('run')) %>% filter(type == 'Ref', run1_EB9223_GCCAAT > 10) %>% pull(run1_EB9215_GATCAG) %>% sd
# [1] 54.69438
# 
# > pltMPRA %>% filter(rs == 'ALAS21') %>% spread(sample, count) %>% dplyr::select(type, contains('run')) %>% filter(type == 'Ref', run1_EB9223_GCCAAT > 10) %>% pull(run1_EB9215_GATCAG) %>% var
# [1] 2991.475
# > ?dnbinom
# > (2991.475 - 86.59375)**-1
# [1] 0.0003442482
# > (2991.475 - 86.59375)**-1 / 86.59375**2
# [1] 4.590907e-08
# > (2991.475 - 86.59375)**-1 * 86.59375**2
# [1] 2.581337

pltMPRA %>% filter(rs == 'ALAS21') %>% spread(sample, count) %>% dplyr::select(type, contains('run')) %>% filter(type == 'Ref', run1_EB9223_GCCAAT > 10) %>% ggplot(aes(run1_EB9215_GATCAG)) + geom_rug() + stat_function(fun = 'dnbinom', args = list(size = 2.581337, mu = 86.59375), xlim = c(0, 250), n = 251)