# Honestly I should probably just do this in one whole tidyverse shot starting
# from an input fasta set + barcode descriptions. Would be cleaner to do that
# way. I suppose.

library(tidyverse)
library(rstan)
select = dplyr::select
filter = dplyr::filter

evalPosterior = function(constrDat, optimalBandwidths, byConstructSummaryStats){
  # Evaluates the posterior for one variant, dependent on the the information
  # empirically obtained from other constructs in the assay. This is done
  # through conditional density estimation on the transcriptional shifts
  # observed in other variants, weighted according to their variance and
  # similarity in annotation space.
  
  lik.grid = seq(-10, 10, by = .1) %>% #TODO trim this grid to a small number of values with more at higher gradient regions. Like creating 3D meshes for video games!
    expand.grid(.,.) %>% #LOOK I MADE AN OWL
    as.tbl %>% 
    mutate(lik.val = map2_dbl(Var1, 
                              Var2, 
                              log.lik.fun, 
                              dat = constrDat, 
                              sig = 1.2 * (byConstructSummaryStats %>% filter(construct == constrDat$construct[1]) %>% .$pooledSD)))
  
  lik.grid %>% 
    ggplot(aes(Var1, Var2)) + 
    geom_raster(aes(fill = lik.val), alpha = .667) + 
    geom_contour(aes(z = lik.val), color = 'grey90') +
    xlab('mean mutant activity') + 
    ylab('mean reference activity') + 
    ggtitle(paste0('log-Likelihood function for construct chr', 
                   constrDat$chr[1],
                   ' position ',
                   constrDat$pos[1],
                   ' ',
                   constrDat$ref[1], 
                   ' --> ',
                   constrDat$alt[1])) +
    geom_abline(slope = 1, intercept = 0, color = 'grey50') + 
    geom_segment(arrow = arrow(length = unit(.16, 'cm')), aes(x = 0, y = 0, xend = .5, yend = -.5)) + 
    annotate(geom = 'text', label = 'TS > 0', x = .25, y = -.75, size = 3)
  
  
  
}

getConstructSummaryStats = function(constrDat){
  # return observed TS, means, and pooled sigma for a given construct
  res = constrDat %>% 
    summarize(chr = constrDat$chr[1],
              pos = constrDat$pos[1],
              ref = constrDat$ref[1],
              alt = constrDat$alt[1],
              muRef = mean(qnact[type == 'Ref']),
              muMut = mean(qnact[type == 'Mut']),
              observedTS = muMut - muRef,
              typeTable = list(table(type)), 
              nRef = map_int(typeTable, ~.['Ref']),
              nMut = map_int(typeTable, ~.['Mut']),
              pooledSD = ((nRef - 1)*sd(qnact[type == 'Ref'])^2 + (nMut - 1)*sd(qnact[type == 'Mut'])^2) / (nRef + nMut - 2))
  return(res)
}

optimizeBandwidths = function(byConstructSummaryStats, annotationSources){
  # From the observed data and annotation sources
  # I'm undecided whether this function should look up the annotations itself or load them from a preloaded source somewhere
}

bayesianMPRA = function(MPRA.qnactivity){
  # qnMPRAdat = a data_frame containing all of the data from a MPRA assay (like MPRA.qnactivity)
  #    * Needs to be in a tidy format - one row per construct
  
  #First take out constructs that don't have both reference and mutant alleles
  goodConstructs = MPRA.qnactivity %>% 
    group_by(construct) %>% 
    summarise(typeTable = list(table(type %>% as.vector)), 
              nTypes = map_int(typeTable, length)) %>% 
    filter(nTypes > 1)
  
  MPRA.qnactivity %<>% filter(construct %in% goodConstructs$construct)
  
  # Then compute summary statistics for each construct that are used downstream (observed TS, variance, etc)
  byConstructSummaryStats = MPRA.qnactivity %>% 
    group_by(construct) %>% 
    by_slice(getConstructSummaryStats, .collate = 'rows') %T>% 
    save(file = paste0('outputs/byConstructSummaryStats.RData'))
  
  
  optimalBandwidths = optimizeBandwidths(byConstructSummaryStats, annotationSources)
  
  res = MPRA.qnactivity %>% 
    group_by(construct) %>% 
    evalPosterior(optimalBandwidths, byConstructSummaryStats)
  return(res)
}

### Likelihood function ####
log.lik.fun = function(muMut, muRef, sig, dat){
  
  #dat - data frame containing a data on a single construct (i.e. MPRA.qnactivity in data/gatheredMPRA.RData)
  #returns the mean cross-validated likelihood of the (h1,h2) bandwidth pair
  
  library(purrr)
  library(dplyr)
  return(dat %>% 
           summarise(logLikRef = dnorm(qnact[type == 'Ref'] - muRef, 
                                       mean = 0,
                                       sd = sig,
                                       log = TRUE) %>% sum, 
                     logLikMut = dnorm(qnact[type == 'Mut'] - muMut,
                                       mean = 0,
                                       sd = sig,
                                       log = TRUE) %>% sum,
                     logLikTot = logLikRef + logLikMut) %>%
           .$logLikTot )
  
}

posterior = 'Insert STAN code here'


### Junk code - kept for reference ####
# getConstructSummaryStatsWithData = function(constrDat){
#   # return observed TS, means, and pooled sigma for a given construct
#   # This returns the data iself as a list column. This doesn't work however as
#   # for some reason nested tibble columns consume huge amounts of memory. With
#   # this data anyway.
#   smallDat = constrDat %>% select(barcode:qnact)
#   res = constrDat %>% 
#     summarize(constrDat = list(smallDat),
#               ref = constrDat$ref[1],
#               alt = constrDat$alt[1],
#               muRef = mean(qnact[type == 'Ref']),
#               muMut = mean(qnact[type == 'Mut']),
#               observedTS = muMut - muRef,
#               typeTable = list(table(type)), 
#               nRef = map_int(typeTable, ~.['Ref']),
#               nMut = map_int(typeTable, ~.['Mut']),
#               pooledSD = ((nRef - 1)*sd(qnact[type == 'Ref'])^2 + (nMut - 1)*sd(qnact[type == 'Mut'])^2) / (nRef + nMut - 2))
#   return(res)
# }
