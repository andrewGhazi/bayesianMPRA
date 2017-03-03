#Now let's try to derive a conditional density estimate prior using 11 histone marks
library(tidyverse)
load('data/varInfoWithHistoneMarkAnnotations.RData')
names(varInfo)[12] = 'transcriptionalShift'

#### Optimize the bandwidths h1 & hN ####
intLeaveOneOutDens = function(y, x, h1, hN){
  
}

AMISE = function(y, x, h1, hN){
  library(purrr) #for map_dbl()
  N = length(y)
  firstTerm = 1/(length(y)) * sum(map_dbl(1:N, intLeaveOneOutDens))
  secondTerm = 
  
  return(firstTerm - secondTerm)
}

#### For one variant ####
prior.fun = function(h1, hN, annotations, TS){
  #sum normal distributions centered at TS values of other variants in the data weighted by optimized h1, hN, and their own standard deviations
  
}