library(tidyverse)

dir = "/mnt/labhome/andrew/MPRA/paper_data/"

UMPRA = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                   delim = "\t",
                   col_names = T,
                   col_types = cols(chr = "c"))

depthNormalize = function(sampleCounts){
  sampleCounts*1e6/sum(sampleCounts)
}

barcode_activities = UMPRA %>% 
  mutate_at(vars(contains('K562')), depthNormalize) %>% 
  mutate(dnaMean = (K562_minP_DNA1 + K562_minP_DNA2)/2) %>% 
  filter(dnaMean > .13) %>% 
  dplyr::select(-K562_minP_DNA1, -K562_minP_DNA2) %>% 
  dplyr::select(-matches('GATA')) %>% 
  gather(sample, depthAdjCount, K562_CTRL_minP_RNA1:K562_CTRL_minP_RNA6) %>% 
  mutate(activity = log(depthAdjCount / dnaMean)) 

construct_data = barcode_activities %>% 
  group_by(construct) %>% 
  nest()

construct_tests = construct_data %>% 
  filter(map_lgl(data, ~sum(.x$type == 'Ref') > 0 & sum(.x$type == 'Mut') > 0)) %>%
  mutate(u_test = map(data, ~wilcox.test(x = .x$activity[.x$type == 'Ref' & is.finite(.x$activity)], y = .x$activity[.x$type == 'Mut' & is.finite(.x$activity)])),
         u_p = map_dbl(u_test, ~.x$p.value), 
         u_q = p.adjust(u_p, method = 'fdr'),
         t_test = map(data, ~t.test(x = .x$activity[.x$type == 'Ref' & is.finite(.x$activity)], y = .x$activity[.x$type == 'Mut' & is.finite(.x$activity)])),
         t_p = map_dbl(t_test, ~.x$p.value),
         t_q = p.adjust(t_p, method = 'fdr'),
         ts_estimate = map_dbl(t_test, ~.x$estimate['mean of y'] - .x$estimate['mean of x']))

t_functional = construct_tests %>% 
  filter(t_q < .01)

u_functional = construct_tests %>% 
  filter(u_q < .01)

# The u and t-test mostly agree (58 shared out of 67 by t and 66 by u)
dplyr::intersect(t_functional %>% select(-data,  -u_test, -t_test), 
                 u_functional %>% select(-data, -u_test, -t_test))

t_functional = construct_tests %>% 
  filter(t_q < .01)

u_functional = construct_tests %>% 
  filter(u_q < .01)

# The u and t-test mostly agree (58 shared out of 67 by t and 66 by u)
test_consistent = dplyr::intersect(t_functional %>% select(-data,  -u_test, -t_test), 
                 u_functional %>% select(-data, -u_test, -t_test))

test_consistent %>%
  ggplot(aes(ts_estimate)) +
  geom_histogram() + 
  geom_rug() +
  scale_x_continuous(limits = c(-3,3))