n = 400
jitterWidth = .2

dat = data_frame(type = c(rep('Ref', 40), rep('Mut', 40), rep('Ref', n), rep('Mut', n)), 
                 Activity = c(rnorm(40, 0, 1.6), rnorm(40, 1.2, 1.6), rnorm(n, 0, 1.6), rnorm(n, 1.2, 1.6)),
                 source = c(rep('Original Variant', 2*40), rep('Similar Variants', 2*n)),
                 xpos = c(rep(0, 40) + runif(40, min = -jitterWidth, max = jitterWidth),
                          rep(1, 40) + runif(40, min = -jitterWidth, max = jitterWidth),
                          rep(0, n) + runif(n, min = -jitterWidth, max = jitterWidth),
                          rep(1, n) + runif(n, min = -jitterWidth, max = jitterWidth)))

dat %>% 
  filter(source == 'Original Variant') %>% 
  ggplot(aes(xpos, Activity, alpha = source)) + 
  geom_point() +
  ylim(range(dat$Activity)) +
  scale_alpha_discrete(range = c(1, .2)) +
  theme_bw() +
  xlab('Type') +
  scale_x_continuous(breaks = c(0,1),
                     limits = c(-.5, 1.5),
                     labels = c('Mut', 'Ref')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20))
#export at 800x550

dat %>% 
  ggplot(aes(xpos, Activity, alpha = source)) + 
  geom_point() +
  ylim(range(dat$Activity)) +
  scale_alpha_discrete(range = c(1, .2)) +
  theme_bw() +
  xlab('Type') +
  scale_x_continuous(breaks = c(0,1),
                     limits = c(-.5, 1.5),
                     labels = c('Mut', 'Ref')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20))

MPRA.qnactivity %>% 
  filter(construct == '10 46044857 1/3') %>% 
  ggplot(aes(type, qnact)) + 
  geom_jitter(height = 0, width = .25) +
  ggtitle('Failure Case: chr10:46044857 1/3 G --> A') +
  theme(text = element_text(size = 20)) +
  ylab('Quantile Normalized Activity')
