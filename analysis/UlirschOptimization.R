#Okay, let's try to do the bandwidth optimization on the Ulirsch data

library(dplyr)
library(magrittr)
library(parallel)
library(ggplot2)

load('~/MPRA/data/varInfo.RData')

#Oh God why haven't I saved this when I did it before
# load('~/MPRA/data/gatheredMPRA.RData')
# 
# poolsig <-function(dat){
#   sds = summarise(group_by(dat, type), sd = sd(qnact))
#   
#   smut = sds$sd[sds$type == 'Mut']
#   sref = sds$sd[sds$type == 'Ref']
#   nref = sum(dat$type == 'Ref')
#   nmut = sum(dat$type == 'Mut')
#   
#   return(sqrt(((nref - 1)*sref**2 + (nmut - 1)*smut**2) / (nref + nmut - 2)))
# }
# 
# varInfo$pooledSigma = mclapply(varInfo$construct, function(x){
#   mpradat = filter(MPRA.qnactivity, construct == x)
#   return(poolsig(mpradat))
# }, mc.cores = 21) %>% unlist

# Optimize bandwidths

# Run model
