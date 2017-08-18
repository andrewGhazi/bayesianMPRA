library(tidyverse)
library(magrittr)
library(rstan)
library(parallel)
library(fitdistrplus)
source('R/mledistModified.R')

### Stan model -----
modelString = "
data{
  int<lower=0> nRefBarcode ; // number of barcodes in ref allele
  int<lower=0> nMutBarcode ;
  int<lower=0> nDNAblocks ; // number of DNA replicates / blocks / samples / transfections
  int<lower=0> nRNAblocks ;
  int<lower=0> refDNAmat[nRefBarcode, nDNAblocks] ; // MPRA count matrix. Rows = barcodes, columns = samples/blocks
  int<lower=0> refRNAmat[nRefBarcode, nRNAblocks] ;
  int<lower=0> mutDNAmat[nMutBarcode, nDNAblocks] ;
  int<lower=0> mutRNAmat[nMutBarcode, nRNAblocks] ;
  real<lower=0> muRefRNAHyperParams[2] ; // gamma hyper-parameters on negative binomial parameters
  real<lower=0> phiRefRNAHyperParams[2] ;
  real<lower=0> muMutRNAHyperParams[2] ;
  real<lower=0> phiMutRNAHyperParams[2] ;
  real<lower=0> muDNAHyperParams[2] ;
  real<lower=0> phiDNAHyperParams[2] ;
}
parameters{
  real<lower=0> muRefDNA[nDNAblocks] ; //mean parameters for each block in each allele for each nucleic acid
  real<lower=0> muRefRNA[nRNAblocks] ;
  real<lower=0> muMutDNA[nDNAblocks] ;
  real<lower=0> muMutRNA[nRNAblocks] ;
  real<lower=0> phiRefDNA[nDNAblocks] ; // size parameters
  real<lower=0> phiRefRNA[nRNAblocks] ;
  real<lower=0> phiMutDNA[nDNAblocks] ;
  real<lower=0> phiMutRNA[nRNAblocks] ;
}
model{
  
  // negative binomial parameters come from gamma hyper-priors
  muRefDNA ~ gamma(muDNAHyperParams[1], muDNAHyperParams[2]) ; 
  muMutDNA ~ gamma(muDNAHyperParams[1], muDNAHyperParams[2]) ;
  phiRefDNA ~ gamma(phiDNAHyperParams[1], phiDNAHyperParams[2]) ;
  phiMutDNA ~ gamma(phiDNAHyperParams[1], phiDNAHyperParams[2]) ;
  
  muRefRNA ~ gamma(muRefRNAHyperParams[1], muRefRNAHyperParams[2]) ;
  muMutRNA ~ gamma(muMutRNAHyperParams[1], muMutRNAHyperParams[2]) ;
  
  phiRefRNA ~ gamma(phiRefRNAHyperParams[1], phiRefRNAHyperParams[2]) ;
  phiMutRNA ~ gamma(phiMutRNAHyperParams[1], phiMutRNAHyperParams[2]) ;

  // count data comes from the specified negative binomial
  for (i in 1:nDNAblocks){ 
    refDNAmat[,i] ~ neg_binomial_2(muRefDNA[i], phiRefDNA[i]) ;
    mutDNAmat[,i] ~ neg_binomial_2(muMutDNA[i], phiMutDNA[i]) ;
  }

  for (i in 1:nRNAblocks){
    refRNAmat[,i] ~ neg_binomial_2(muRefRNA[i], phiRefRNA[i]) ;
    mutRNAmat[,i] ~ neg_binomial_2(muMutRNA[i], phiMutRNA[i]) ;
  }
}
"

model = stan_model(model_code = modelString)

set.seed(1280)

### By variant prior estimation ----
load('data/varInfoWithHistoneMarkAnnotations.RData') 
names(varInfo)[12] = 'transcriptionalShift'

varInfo %<>% dplyr::select(construct, chr:alt, gkmerDelta:transcriptionalShift)

# varInfo %>% select(contains('Broad'), contains('Sydh'), eigen:DeepSeaDnaase)
preds = varInfo %>% 
  dplyr::select(DeepSeaDnaase, gkmerDelta) %>%
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

distMat = generateDistMat(preds)

findWeights = function(i, distMat, minDistKernel, minNumContributing = 30, increaseFold = 1.333){
  # for the ith variant, produce a vector of weightings such that at minimum minNumContributing
  # variants meaningfully contribute (i.e. cumsum(sortedWeights)[30] <= .99)
  # arguments after the 2nd are heuristics that may need tuning
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix (precomputed for speed)
  distKernel = minDistKernel
  
  rawWeights = dnorm(distMat[i,-i], sd = distKernel)
  scaledWeights = rawWeights / sum(rawWeights)
  sorted = sort(scaledWeights, decreasing = TRUE)
  
  notEnoughContributing = cumsum(sorted[1:minNumContributing])[minNumContributing] > .99
  #allZero = all(rawWeights$x == 0)
  
  # if there aren't more than minNumContributing variants providing meaningful contribution to the prior
  if (is.na(notEnoughContributing) || notEnoughContributing) { 
    
    # iteratively increase the kernel bandwith until they do
    while (is.na(notEnoughContributing) || notEnoughContributing) {
      distKernel = distKernel * increaseFold
      rawWeights = dnorm(distMat[i,-i], sd = distKernel)
      scaledWeights = rawWeights / sum(rawWeights)
      sorted = sort(scaledWeights, decreasing = TRUE)
      
      notEnoughContributing = cumsum(sorted[1:minNumContributing])[minNumContributing] > .99
    }
  }
  
  scaledWeights
}

# Initialize the kernel at some small value based on the typical distances in the input distance matrix
minDistKernel = distMat[upper.tri(distMat)] %>% 
  unlist() %>% 
  sort() %>% #sort all observed distances
  .[. > 0] %>% 
  quantile(.001) # pick the .1th quantile. The only variants that will use this kernel will be in very densely populated regions of predictor space

varInfo %<>% 
  mutate(weights = map(1:nrow(.), ~findWeights(.x, distMat, minDistKernel)))

### read in counts -----------
dir = "/mnt/labhome/andrew/MPRA/paper_data/"
ulirschCounts = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                           delim = "\t",
                           col_names = T,
                           col_types = cols(chr = "c")) %>% 
  dplyr::select(matches('construct|CTRL|DNA|type')) %>% 
  group_by(construct) %>% 
  nest(.key = countData)

varInfo %<>% left_join(ulirschCounts, by = 'construct')

### Fit NegBin estimates to each variant -----
safelyFitNegBin = function(countVec){
  #This doesn't quite suppress all error messages but it does output the right things
  #an error can still be printed when there's very low variability of low counts (e.g. c(1, rep(0, 14)))
  safely(fitdist, 
         otherwise = list(estimate = purrr::set_names(rep(NA, 2), nm = c('mu', 'size'))), 
         quiet = TRUE)(countVec, 'nbinom')
}

estTransfectionParameters = function(countDat){
  # countDat - a tibble with a type (ref/mut) column and columns of observed MPRA counts in given transfections
  # uses fitdistrplus::fitdist because MASS::fitdistr was cracking wise at me
  # using a modified version of fitdistrplus::mledist because the regular version can't use non-integer weights
  
  countDat %>% 
    gather(block, count, -type) %>% 
    group_by(type, block) %>% 
    summarise(MLEnegBin = list(safelyFitNegBin(count))) %>% 
    ungroup %>% 
    mutate(muEst = map_dbl(MLEnegBin, ~.x$result$estimate['mu']),
           sizeEst = map_dbl(MLEnegBin, ~.x$result$estimate['size'])) %>% 
    dplyr::select(-MLEnegBin)
}


ncores = 20

varInfo %<>% 
  mutate(negBinParams = mclapply(countData, estTransfectionParameters, mc.cores = ncores))

paramDF = varInfo %>% dplyr::select(construct, negBinParams)

### Fit weighted gamma hyperprior on negBin parameters to each variant ------- TODO fix this code
fitMuGamma = function(weights, muEstimates){
  fnToMin = function(paramVec){-sum(weights * dgamma(muEstimates, 
                                                     shape = paramVec[1], 
                                                     rate = paramVec[2], 
                                                     log = TRUE))}
  
  meanEst = mean(muEstimates)
  varEst = var(muEstimates)
  
  initialGuess = c(meanEst**2 / varEst, meanEst/varEst)
  
  optimRes = optim(initialGuess, 
                   fnToMin, 
                   lower = c(1e-12), 
                   control = list(ndeps = c(1e-4, 1e-5)))
  
  # The parscale control option is needed to scale the optimization proposals.
  # If the rate suggestions are too huge then fnToMin throws out Inf which
  # breaks the optimizer
  
  if (optimRes$convergence != 0) {
    stop(paste0('problems with gamma fitting, convergence code: ', optimRes$convergence))
  }
  
  optimRes$par %>% 
    set_names(c('shape', 'rate'))
}

fitSizeGamma = function(weights, sizeEstimates){
  # different ndeps, no lower bound so the optimizer works
  fnToMin = function(paramVec){-sum(weights * dgamma(sizeEstimates, 
                                                     shape = paramVec[1], 
                                                     rate = paramVec[2], 
                                                     log = TRUE))}
  meanEst = mean(sizeEstimates)
  varEst = var(sizeEstimates)
  
  initialGuess = c(meanEst**2 / varEst, meanEst/varEst)
  
  optimRes = optim(initialGuess, 
                   fnToMin, 
                   control = list(ndeps = c(1e-4, 1e-4))) 
  
  if (optimRes$convergence != 0) {
    stop(paste0('problems with gamma fitting, convergence code: ', optimRes$convergence))
  }
  
  optimRes$par %>% 
    set_names(c('shape', 'rate'))
}

plusOrHomebrewMu = function(weights, muEstimates, initialMuGuess){
  # Try to fit the gamma with fitdistrplus. If that doesn't work, try the
  # homebrew optimizer that calculates its own type-block initial guess for the
  # individual type-block
  
  res = try(fitdistMod(muEstimates, 
                       'gamma', 
                       weights = weights, 
                       start = initialMuGuess,
                       control = list(ndeps = c(1e-4, 1e-5)),
                       lower = 1e-12)$estimate,
            silent = TRUE)
  
  if (class(res) == 'try-error') {
    fitMuGamma(weights, muEstimates)
  } else{
    res
  }
}

plusOrHomebrewSize = function(weights, sizeEstimates, initialSizeGuess){
  #same as above just different ndeps
  
  res = try(fitdistMod(sizeEstimates, 
                       'gamma', 
                       weights = weights, 
                       start = initialSizeGuess,
                       control = list(ndeps = c(1e-4, 1e-4)),
                       lower = 1e-12)$estimate,
            silent = TRUE)
  
  if (class(res) == 'try-error') {
    fitSizeGamma(weights, sizeEstimates)
  } else{
    res
  }
}

fitGammaHyperPriors = function(constructNum){ 
  constr = varInfo[constructNum,]
  others = varInfo[-constructNum,]
  
  
  dataForEstimates = others %>% 
    mutate(weight = constr$weights[[1]]) %>% 
    dplyr::select(construct, negBinParams, weight) %>% 
    unnest %>% 
    na.omit %>% 
    filter(weight > 1e-4*sort(constr$weights[[1]], decreasing = TRUE)[30]) # This is necessary for speed)
  
  muDat = dataForEstimates$muEst
  initialMuGuess = list(shape = mean(muDat)**2 / var(muDat), 
                        rate = mean(muDat) / var(muDat))
  
  sizeDat = dataForEstimates$sizeEst[dataForEstimates$sizeEst < quantile(dataForEstimates$sizeEst, .99)]
  initialSizeGuess = list(shape = mean(sizeDat)**2 / var(sizeDat), 
                        rate = mean(sizeDat) / var(sizeDat))
  
  dataForEstimates %>% 
    group_by(type, block) %>% 
    summarise(muGammaHyperPriors = list(plusOrHomebrewMu(weight, 
                                                         muEst, 
                                                         initialMuGuess)),
              sizeGammaHyperPriors = list(plusOrHomebrewSize(weight[sizeEst < 1e4],
                                                             sizeEst[sizeEst < 1e4], 
                                                             initialSizeGuess))) %>% 
    ungroup
  
  # Estimates of the size parameter can be unstable (because variance = mu + mu^2
  # / size so size --> Inf as var --> mu) so we cut out those that are above the
  # 99th quantile. These are variants that are essentially poisson in their
  # counts so we're slightly biasing our result to show HIGHER variance.
}

varInfo %<>% 
  mutate(gammaParams = mclapply(1:n(), fitGammaHyperPriors, mc.cores = 20))

# system.time(varInfo <- varInfo %>%
#               mutate(gammaParams = mclapply(1:n(), fitGammaHyperPriors, mc.cores = 20)))

save(varInfo, file = '~/bayesianMPRA/outputs/varInfoWithNegBinAndGammaParams.RData')

varInfo[varInfo$kNN[[1]],]$negBinParams %>% 
  purrr::reduce(bind_rows) %>% 
  group_by(type, block) %>% 
  summarise(muHyperParams = list(safelyFitGamma(muEst)),
            sizeHyperParams = list(safelyFitGamma(sizeEst))) %>% 
  ungroup

## TODO: keep adapting over ulirschNegBinPrior ------
## 1. fit neg binomials
## 2. fit WEIGHTED gamma hyperprior to RNA counts
## 3. Fit marginal gamma hyperprior on DNA counts
## 4. adapt model code
## 5. run sampling
## 6. money
