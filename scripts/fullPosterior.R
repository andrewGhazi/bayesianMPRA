# Honestly I should probably just do this in one whole tidyverse shot starting
# from an input fasta set + barcode descriptions. Would be cleaner to do that
# way. I suppose.

library(tidyverse)
library(rstan)
select = dplyr::select
filter = dplyr::filter

evalPosterior = function(constrDat, otherConstrDat){
  # Evaluates the posterior for one variant, dependent on the the information
  # empirically obtained from other constructs in the assay
  
  
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

optimizeBandwidths = function(byConstructSummaryStats){
  # From the observed data and annotation sources
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
  
  
  
  optimalBandwidths = optimizeBandwidths(byConstructSummaryStats, annotationD)
  
  res = MPRA.qnactivity %>% 
    group_by(construct) %>% 
    evalPosterior
  return(res)
}

### Likelihood function ####
log.lik.fun = function(dat, h1, hN, TS){
  
  
  #dat - data frame containing a data on a single construct (i.e. MPRA.qnactivity in data/gatheredMPRA.RData)
  # h1 = bandwidth value 1 - along transcriptional shift
  # hN = a vector of the N bandwidths for each annotation (i.e. one for each column of annotations)
  #returns the mean cross-validated likelihood of the (h1,h2) bandwidth pair
  
  # input checks
  if (nrow(annotations) != length(y)) {stop('There aren\'t as many annotation rows as there are transcriptional shifts')}
  if (ncol(annotations) != length(hN)) {stop('There aren\'t as many annotation bandwidths as there are in the annotation matrix')}
  if (names(y) != c('transcriptionalShift', 'refActStdDev')) {
    stop('Please input the transcriptional shift and reference construct activity standard deviations as a data_frame with columns titled \'transcriptionalShift\' and \'refActStdDev\'')
  }
  
  library(purrr)
  library(dplyr)
  
  getVarLogLik = function(varDat, TS){
    
    return(varDat %>% 
             summarise(logLik = dnorm(qnact[type == 'Mut'] - TS, 
                                      mean = mean(qnact[type == 'Ref']), 
                                      sd = sd(qnact[type == 'Ref']) , 
                                      log = TRUE) %>% sum()) %>% .$logLik)
  }
  
  dat %>% group_by(construct) %>% by_slice(getVarLogLik, TS = 1, .collate = 'rows') %>% .$.out #don't need to sum because we're evaluating each variant separately
  
}

posterior = 'Insert STAN code here'
