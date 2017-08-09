library(tidyverse)
library(magrittr)
library(rstan)

### Stan model -----
modelString = "
data{
int<lower=0> nRefBarcode ;
int<lower=0> nMutBarcode ;
int<lower=0> nDNAblocks ;
int<lower=0> nRNAblocks ;
int<lower=0> refDNAmat[nRefBarcode, nDNAblocks] ;
int<lower=0> refRNAmat[nRefBarcode, nRNAblocks] ;
int<lower=0> mutDNAmat[nMutBarcode, nDNAblocks] ;
int<lower=0> mutRNAmat[nMutBarcode, nRNAblocks] ;
real<lower=0> muHyperParams[2] ;
real<lower=0> phiHyperParams[2] ;
}
parameters{
real<lower=0> muRefDNA[nDNAblocks] ;
real<lower=0> muRefRNA[nRNAblocks] ;
real<lower=0> muMutDNA[nDNAblocks] ;
real<lower=0> muMutRNA[nRNAblocks] ;
real<lower=0> phiRefDNA[nDNAblocks] ;
real<lower=0> phiRefRNA[nRNAblocks] ;
real<lower=0> phiMutDNA[nDNAblocks] ;
real<lower=0> phiMutRNA[nRNAblocks] ;
}
model{

muRefDNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;
muRefRNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;
muMutDNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;
muMutRNA ~ gamma(muHyperParams[1], muHyperParams[2]) ;

phiRefDNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
phiRefRNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
phiMutDNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;
phiMutRNA ~ gamma(phiHyperParams[1], phiHyperParams[2]) ;

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

findWeights = function(i, distMat, minNumContributing = 30){
  # for the ith variant, produce a vector of weightings such that at minimum minNumContributing
  # variants meaningfully contribute (i.e. cumsum(sortedWeights)[30] <= .99)
  
  distKernel = .01
  rawWeights = sort(dnorm(distMat[i,-i], sd = distKernel), decreasing = TRUE, index.return = TRUE)
  scaledWeights = rawWeights$x / sum(rawWeights$x)
  
  # if there aren't more than minNumContributing variants providing meaningful contribution to the prior
  if (cumsum(scaledWeights[1:minNumContributing])[minNumContributing] > .99) { #minNumContributing and .99 are both heuristics that may need tuning
    
    # iteratively increase the kernel bandwith until they do
    while (cumsum(scaledWeights[1:minNumContributing])[minNumContributing] > .99) {
      distKernel = distKernel * 1.5
      rawWeights = sort(dnorm(distMat[i,-i], sd = distKernel), decreasing = TRUE, index.return = TRUE)
      scaledWeights = rawWeights$x / sum(rawWeights$x)
    }
  }
  
  scaledWeights
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
varInfo %<>% 
  mutate(kNN = map(1:nrow(.), ~findKNN(k, .x, distMat)),
         numNeigh = map_int(kNN, length))

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