library(tidyverse)
library(stringr)
library(magrittr)

load('~/bayesianMPRA/outputs/varFuns.RData')

HDIofMCMC = function(sampleVec, credMass=0.95) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
} # From Kruschke

loadAndDo = function(construct) {
  #loads once and does everything necessary
  fileName = construct %>% str_replace_all(' ', '_') %>% str_replace('/', '-') %>% str_c('.RData')
  load(paste0('~/bayesianMPRA/outputs/30NNmcmcOutputs/', fileName))
  
  data_frame(construct = construct) %>% 
    transmute(
              ts95HDI = HDIofMCMC(res$batch[,2] - res$batch[,1], credMass = .95) %>% list,
              ts99HDI = HDIofMCMC(res$batch[,2] - res$batch[,1], credMass = .99) %>% list,
              probablyFunctional = map_lgl(ts95HDI, ~!(.x[1] < 0 && .x[2] > 0)),
              definitelyFunctional = map_lgl(ts99HDI, ~!(.x[1] < 0 && .x[2] > 0)),
              meanRefMu = mean(res$batch[,1]),
              meanMutMu = mean(res$batch[,2]),
              meanSig = mean(res$batch[,3]),
              adNormTestPValue = goftest::ad.test(res$batch[,2] - res$batch[,1], null = 'pnorm')$p.value) # This doesn't work, needs to be fixed.

}

varFuns %<>% rowwise %>% bind_cols(., do(., loadAndDo(.$construct))) %>% ungroup

save(varFuns, file = '~/bayesianMPRA/outputs/varFuns.RData')
