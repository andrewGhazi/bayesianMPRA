library(tidyverse)
library(magrittr)
library(stringr)
library(mcmc)
library(parallel)

select = dplyr::select

load('data/varInfoWithHistoneMarkAnnotations.RData') 
load('data/gatheredMPRA.RData')
names(varInfo)[12] = 'transcriptionalShift'

# varInfo %>% select(contains('Broad'), contains('Sydh'), eigen:DeepSeaDnaase)
preds = varInfo %>% select(DeepSeaDnaase, gkmerDelta) %>%
  map_df(~scale(.x)[,1]) #effing scale only outputs matrices

generateDistMat = function(predictors) {
  #predictors is a n x d data frame of predictors
  
  predictors %>% 
    as.data.frame() %>% 
    as.matrix %>% 
    dist %>% 
    as.matrix %>% 
    log1p
  
}

findKNN = function(k, i, distMat){
  # For the ith variant, return the indices of the k nearest neighbors in the dist mat
  
  res = sort(distMat[i,], index.return = TRUE)
  
  # If there are a lot that have exactly the same predictors, just return the indices of all of them
  if (all(res$x[1:(k+1)] == 0)) {# k + 1 because it will always return the 0 for the ith item
    return(res$ix[res$x == 0])
  } 
  
  res$ix = res$ix[!(res$ix %in% i)]
  res$ix[1:k]
}

k = 30 
distMat = generateDistMat(preds)
#save(distMat, file = '~/bayesianMPRA/outputs/distMat.RData')
# load('~/bayesianMPRA/outputs/distMat.RData')

varInfo %<>% 
  mutate(kNN = map(1:nrow(.), ~findKNN(k, .x, distMat)),
         numNeigh = map_int(kNN, length))

likFunFromData = function(refAct, mutAct){
  # Returns a log-likelihood function from two vectors of activities
  likfun =  function(paramVec) {
    # Formerly took three arguments: refMu, mutMu, sig. But mcmc::metrop() needs it to take a vector, so it's changed accordingly
    sum(c(dnorm(refAct, paramVec[1], paramVec[3], log = TRUE), 
          dnorm(mutAct, paramVec[2], paramVec[3], log = TRUE)))
  }
  
  return(likfun)  
}

priorFromLikFuns = function(likFuns, paramVec){ #the prior will be the GEOMETRIC mean of the likelihood functions of the kNN
  priorfun = function(paramVec) {
    # Formerly took three arguments: refMu, mutMu, sig. But mcmc::metrop() needs it to take a vector, so it's changed accordingly
    mean(map_dbl(likFuns, ~.x(paramVec)))
  }
  
  return(priorfun)
}

getPriorFun = function(kNN){
  priorFromLikFuns(varFuns$varLikFun[kNN], paramVec)
}

varFuns = MPRA.qnactivity %>% 
  filter(grepl('CTRL', Block)) %>% 
  group_by(construct) %>% 
  summarise(refDat = list(qnact[type == 'Ref']),
            mutDat = list(qnact[type == 'Mut']),
            varLikFun = map2(refDat, mutDat, ~likFunFromData(.x, .y))) %>% 
  select(-refDat, -mutDat) %>% 
  right_join(varInfo, by = 'construct')

varFuns %<>% mutate(priorFun = lapply(varFuns$kNN, getPriorFun)) #gasp, lapply worked where map didn't

getPostFun = function(likFun, priorFun) {
  function(paramVec){
    if(paramVec[3] < 0) return(-Inf)
    likFun(paramVec) + priorFun(paramVec)
  }
}
varFuns %<>% mutate(postFun = lapply(1:nrow(varFuns), function(i){getPostFun(varFuns$varLikFun[[i]],
                                                                             varFuns$priorFun[[i]])}))

# meanDiffs = MPRA.qnactivity %>% 
#   group_by(construct) %>% 
#   summarise(TS = mean(qnact[type == 'Mut']) - mean(qnact[type == 'Ref']),
#             pooledSD = sqrt(((sum(type == 'Mut') - 1)*sd(qnact[type == 'Mut'])**2 + sum(type == 'Ref') - 1)*sd(qnact[type == 'Ref'])**2 / (length(qnact) - 2)))
# 
# varFuns %<>% left_join(meanDiffs, by = 'construct')
# save(varFuns, file = '~/bayesianMPRA/outputs/varFuns.RData')

# system.time({tmp = metrop(varFuns$postFun[[981]], initial = c(0,0,1), nbatch = 5000, blen = 1, nspac = 5,
#              scale = .12)}) # ~28s

runUlirschMCMC = mclapply(seq_along(varFuns$construct),
                          function(i){
                            strt = Sys.time()
                            burn = metrop(varFuns$postFun[[i]],
                                          initial = c(0,0,1),
                                          nbatch = 100,
                                          blen = 1,
                                          scale = .12,
                                          nspac = 5)
                            res = metrop(burn,
                                         nbatch = 5000,
                                         blen = 1,
                                         scale = .12,
                                         nspac = 5)
                            save(res, file = paste0('~/bayesianMPRA/outputs/k30gkmDeepSeaPrior/', varFuns$construct[i] %>% str_replace_all(' ', '_') %>% str_replace('/', '-'), '.RData'))
                            stp = Sys.time() - strt
                            return(stp)
                          }, mc.cores = 20)
