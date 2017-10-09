# So the estimated size parameters for the negative binomial can be huge if the counts are low
varInfo[varInfo$kNN[[1]],]$negBinParams %>% 
  purrr::reduce(bind_rows) %>% 
  group_by(type, block) %>% 
  summarise(muHyperParams = list(safelyFitGamma(muEst)),
            sizeHyperParams = list(safelyFitGamma(sizeEst))) %>% 
  ungroup # errors lead me to find a variant with the following counts which produces a size estimate of 8e5

c(3,0,1,0,2,3,2,0,4,1,2,2,0,2) %>% fitdist(distr = 'nbinom')

# Seriously wtf how does a 5 order of magnitude change of the size parameter have such a negligible effect on the likelihood?
varInfo[varInfo$kNN[[1]],]$countData %>% map(~filter(.x, type == 'Mut')$K562_CTRL_minP_RNA1) %>% .[[21]] %>% dnbinom(size = 8.66e5, mu = 1.5, log = TRUE) %>% sum
varInfo[varInfo$kNN[[1]],]$countData %>% map(~filter(.x, type == 'Mut')$K562_CTRL_minP_RNA1) %>% .[[21]] %>% dnbinom(size = 8.66, mu = 1.5, log = TRUE) %>% sum

# I guess it's mostly because the huge parameter doesn't effect the distribution all that much:
data_frame(x = 0:20, 
           sizeEquals8 = dnbinom(x, size = 8, mu = 1.57),
           sizeEquals8e5 = dnbinom(x, size = 8e5, mu = 1.57)) %>% 
  gather(dist, dens, -x) %>% 
  ggplot(aes(x,dens)) + geom_point(aes(color = dist)) +
  ggtitle('Negative Binomial Distributions for low counts with \nmean = 1.57 and two different size parameters')

# But it does throw off the gamma estimation so we'll need to figure out how to avoid this.
