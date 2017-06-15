# Now let's try to run the hierarchical model but estimate the parameters of the prior from the 30kNN

library(tidyverse)
library(stringr)
library(magrittr)
library(rjags)
library(fitdistrplus)
library(coda)
library(parallel)

# code from kNNprior.R to get the kNN-----------
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
varInfo %<>% 
  mutate(kNN = map(1:nrow(.), ~findKNN(k, .x, distMat)),
         numNeigh = map_int(kNN, length))

# read in counts -----------
dir = "/mnt/labhome/andrew/MPRA/paper_data/"
ulirschCounts = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
           delim = "\t",
           col_names = T,
           col_types = cols(chr = "c")) %>% 
  dplyr::select(matches('construct|CTRL|DNA|type')) %>% 
  group_by(construct) %>% 
  nest(.key = countData)

varInfo %<>% left_join(ulirschCounts, by = 'construct')

# estimate prior parameters

set.seed(1280)

dir = "/mnt/labhome/andrew/MPRA/paper_data/"

exampleData = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                         delim = "\t",
                         col_names = T,
                         col_types = cols(chr = "c")) %>% 
  dplyr::select(matches('construct|CTRL|DNA|type')) %>% 
  filter(construct == '1 155271258 1/3') %>% 
  mutate(bcIndex = 1:nrow(.)) 

neighbors = varInfo %>% 
  filter(construct == '1 155271258 1/3') %>% 
  .$kNN %>% 
  unlist

# NegBin MLE estimates -----------
# For each variant, we fit negative binomial distributions to the counts from each block*allele

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

# Gamma hyper-prior estimation ---------
# Now for each variant, we take the negative binomial parameters of the variant's NEIGHBORS estimated above to define an hyper-prior gamma distribution
# So the 30 nearest neighbor's negative binomials are like 30 samples from the distribution that each block's parameters come from

safelyFitGamma = function(neighborMuEstimates, type, block){
  cat(paste(unique(type), ', ', unique(block), '\n'))
  safely(fitdist)(as.vector(na.omit(neighborMuEstimates)), 'gamma')
}

varInfo[varInfo$kNN[[1]],]$negBinParams %>% 
  purrr::reduce(bind_rows) %>% 
  group_by(type, block) %>% 
  summarise(muHyperParams = list(safelyFitGamma(muEst, type, block)),
            sizeHyperParams = list(safelyFitGamma(sizeEst, type, block))) %>% 
  ungroup


varInfo %>% 
  mutate()

# Now the actual JAGS part --------------

dataList = list(nRefBarcode = sum(exampleData$type == 'Ref'),
                nMutBarcode = sum(exampleData$type == 'Mut'),
                nDNAblocks = colnames(exampleData) %>% grepl('DNA', .) %>% sum,
                nRNAblocks = colnames(exampleData) %>% grepl('RNA', .) %>% sum,
                refDNAmat = exampleData %>% filter(type == 'Ref') %>% dplyr::select(contains('DNA')) %>% as.data.frame %>% as.matrix,
                refRNAmat = exampleData %>% filter(type == 'Ref') %>% dplyr::select(contains('RNA')) %>% as.data.frame %>% as.matrix,
                mutDNAmat = exampleData %>% filter(type == 'Mut') %>% dplyr::select(contains('DNA')) %>% as.data.frame %>% as.matrix,
                mutRNAmat = exampleData %>% filter(type == 'Mut') %>% dplyr::select(contains('RNA')) %>% as.data.frame %>% as.matrix,
                muHyperParams = ,
                sizeHyperParams = )
  
modelStringHier = "
model {
  for ( i in 1:nRefBarcode ) {
    for ( j in 1:nDNAblocks ) {
      refDNAmat[i,j] ~ dnegbin(pDNA[1,j], rDNA[1,j])
    }
    for ( j in 1:nRNAblocks ) {
      refRNAmat[i,j] ~ dnegbin(pRNA[1,j], rRNA[1,j])
    }
  }

  for ( i in 1:nMutBarcode ) {
    for ( j in 1:nDNAblocks ) {
      mutDNAmat[i,j] ~ dnegbin(pDNA[2,j], rDNA[2,j])
    }
    for ( j in 1:nRNAblocks ) {
      mutRNAmat[i,j] ~ dnegbin(pRNA[2,j], rRNA[2,j])
    }
  }

  for ( i in 1:2 ) { #1 = ref, 2 = mut
    for ( j in 1:nDNAblocks) {
      pDNA[i,j] <- rDNA[i,j]/(rDNA[i,j] + mDNA[i,j])
      rDNA[i,j] ~ dgamma(0.01, 0.01)
      mDNA[i,j] ~ dgamma(0.01, 0.01)
      vDNA[i,j] <- rDNA[i,j]*(1-pDNA[i,j])/(pDNA[i,j]*pDNA[i,j])
    }
    for ( j in 1:nRNAblocks) {
      pRNA[i,j] <- rRNA[i,j]/(rRNA[i,j] + mRNA[i,j])
      rRNA[i,j] ~ dgamma(0.01, 0.01)
      mRNA[i,j] ~ dgamma(0.01, 0.01)
      vRNA[i,j] <- rRNA[i,j]*(1-pRNA[i,j])/(pRNA[i,j]*pRNA[i,j])
    }
  }
}
"

#---------------
# adaptive mode adapts sampling parameters to increase efficiency. However this makes it not a markov chain and the samples it provides shouldn't be used
jagsModelHier = jags.model(file = textConnection(modelStringHier),
                           data = dataList,
                           n.chains = 3, 
                           n.adapt = 500) 

#update turns off adaptive mode and provides the actual burn-in
update(jagsModelHier,
       n.iter = 500) 

# and coda.samples uses JAGS to do the actual sampling and output it in a form that's easy to analyze with coda diagnostic tools
codaSamplesHier = coda.samples(jagsModelHier,
                               variable.names = c('mDNA', 'mRNA', 'vDNA', 'vRNA'),
                               n.iter = 3334)

depthScale = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                        delim = "\t",
                        col_names = T,
                        col_types = cols(chr = "c")) %>% 
  dplyr::select(matches('CTRL|DNA')) %>% map_df(sum)

samples = codaSamplesHier %>% 
  map_df(~.x %>% as.matrix %>% as.data.frame %>% as.tbl)

scaleAppropriately = function(dat, paramType, scaleNumber){
  if (paramType == 'm') {
    map_df(dat, ~.x*1e6/scaleNumber) 
  } else {
    map_df(dat, ~.x*1e12/(scaleNumber**2))
  }
}

scaleByDepth = function(colname){
  #Take the column of MCMC parameter samples and scale according to that transfection's depth
  
  acidType = colname %>% str_extract('[DR]NA')
  blockNumber = colname %>% str_extract('[0-9]+]') %>% str_extract('[0-9]')
  paramType = colname %>% str_extract('^[mv]')
  scaleNumber = depthScale %>% 
    dplyr::select(matches(acidType)) %>% 
    dplyr::select(matches(regex(paste0(blockNumber, '$')))) %>% 
    .[[1]] # This will always return a 1x1 df so extracting the value with .[[1]] should be okay
  
  samples %>% 
    .[colname] %>% 
    scaleAppropriately(paramType = paramType, scaleNumber = scaleNumber)
    
}

depthScaledSamples = samples %>% colnames %>% map(scaleByDepth) %>% purrr::reduce(bind_cols)

computeParamBlockMean = function(param){
  depthScaledSamples %>% 
    dplyr::select(contains(param)) %>% 
    rowMeans() %>% 
    as.data.frame %>%
    as.tbl() %>% 
    set_colnames(param)
}


fullPost = c('mDNA[1', 'mDNA[2', 'mRNA[1', 'mRNA[2', 'vDNA[1', 'vDNA[2', 'vRNA[1', 'vRNA[2') %>% #1's = ref, 2's = mut
  map(computeParamBlockMean) %>% 
  purrr::reduce(bind_cols) %>% 
  set_colnames(c('mRefDNA', 'mMutDNA', 'mRefRNA', 'mMutRNA', 'vRefDNA', 'vMutDNA', 'vRefRNA', 'vMutRNA')) %>% 
  mutate(TS = log(mMutRNA/mMutDNA) - log(mRefRNA/mRefDNA))

fullPost %>% 
  ggplot(aes(TS)) + 
  geom_histogram(bins = 50)


