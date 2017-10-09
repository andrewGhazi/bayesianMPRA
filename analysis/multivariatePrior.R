#Now let's try to derive a conditional density estimate prior using 11 histone marks
library(tidyverse)
library(dmvnorm) #I'm using this instead of mvtnorm because the following suggests it is faster: http://stats.stackexchange.com/questions/68476/drawing-from-the-multivariate-students-t-distribution
load('data/varInfoWithHistoneMarkAnnotations.RData')
names(varInfo)[12] = 'transcriptionalShift'

# This creates a group of histograms showing how each predictor is individually distributed. Note that it throws out 40% of the data for the sake of a log scale
# varInfo %>% 
#   select(construct, eigen:DeepSeaDnaase, BroadK562H3k4me1:SydhK562H3k27me3b) %>% 
#   gather(predictor, value, -construct) %>% 
#   ggplot(aes(value)) + 
#   geom_histogram() + 
#   facet_wrap('predictor') + 
#   scale_x_log10()

varInfoAnnotations = varInfo %>% select(construct, PICS:DeepSeaDnaase, BroadK562H3k4me1:SydhK562H3k27me3b)
byConstructSummaryStats %<>% 
  drop_na(-typeTable) %>% 
  left_join(varInfoAnnotations, by = 'construct')

#### Optimize the bandwidths h1 & hN ####
intLeaveOneOutDens = function(y, x, h1, hN){
  
}

leaveOneOutDens = function(y, x, h1, hN, ySD){
  # x needs to alrady be a matrix (for sake of speed)
  N = length(y)
  
  varCov = diag(hN)
  
  for (i in 1:N){
    #use hN & to get weights
    point = y[i]
    rest = y[-i]
    
    pointFPs = x[i,]
    functionalPredictors = x[-i,]
    funPredBasedWeights = dmvnorm(functionalPredictors, mean = pointFPs, sigma = varCov) 
    
    #call density 
    density(y, 
            bw = h1,
            weights = sigAndFunPredBasedWeights)
  }
}

x = byConstructSummaryStats %>% select(BroadK562H3k4me1:SydhK562H3k27me3b) %>% as.data.frame() %>% as.matrix()
AMISE = function(y, x, h1, hN, ySD){
  library(purrr) #for map_dbl()
  N = length(y)
  firstTerm = 1/(length(y)) * sum(map_dbl(1:N, intLeaveOneOutDens))
  secondTerm = 2/(length(y)) * sum(map_dbl(1:N, leaveOneOutDens))
  
  return(firstTerm - secondTerm)
}

#### For one variant ####
prior.fun = function(h1, hN, annotations, TS){
  #sum normal distributions centered at TS values of other variants in the data weighted by optimized h1, hN, and their own standard deviations
  
  #grid search over hN values (crap) for roughly optimized values
}


# mvtnormTime = microbenchmark({funPredBasedWeights = dmvnorm(functionalPredictors, mean = pointFPs, sigma = varCov)}, times = 10000)
# mnormtTime = microbenchmark({funPredBasedWeights = dmnorm(functionalPredictors, mean = pointFPs, varcov = varCov)}, times = 10000)
# mvNormTimes = data_frame(mvtnorm = mvtnormTime$time, mnormt = mnormtTime$time) %>% gather(package, time)
# 
# ggplot(mvNormTimes, aes(package, time)) + geom_violin() + ggtitle('10k trials of evaluating leave-one-out weightings') + scale_y_log10() + ylab('time (ns)')
# ggsave('outputs/multivariateNormWeightTimings.png')
# mvNormTimes %>% group_by(package) %>% summarise(pkgMean = mean(time), pkgMed = median(time), pkgTot = sum(time))
# Takeaway: mvtnorm is just slightly faster on average.
