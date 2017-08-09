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

findWeights = function(i, distMat, minDistKernel, minNumContributing = 30, increaseFold = 1.5){
  # for the ith variant, produce a vector of weightings such that at minimum minNumContributing
  # variants meaningfully contribute (i.e. cumsum(sortedWeights)[30] <= .99)
  # arguments after the 2nd are heuristics that may need tuning
  
  # Initialize the kernel at some small value based on the typical distances in the input distance matrix (precomputed for speed)
  distKernel = minDistKernel
  
  rawWeights = sort(dnorm(distMat[i,-i], sd = distKernel), decreasing = TRUE, index.return = TRUE)
  scaledWeights = rawWeights$x / sum(rawWeights$x)
  
  notEnoughContributing = cumsum(scaledWeights[1:minNumContributing])[minNumContributing] > .99
  #allZero = all(rawWeights$x == 0)
  
  # if there aren't more than minNumContributing variants providing meaningful contribution to the prior
  if (is.na(notEnoughContributing) || notEnoughContributing) { 
    
    # iteratively increase the kernel bandwith until they do
    while (is.na(notEnoughContributing) || notEnoughContributing) {
      distKernel = distKernel * increaseFold
      rawWeights = sort(dnorm(distMat[i,-i], sd = distKernel), decreasing = TRUE, index.return = TRUE)
      scaledWeights = rawWeights$x / sum(rawWeights$x)
      
      notEnoughContributing = cumsum(scaledWeights[1:minNumContributing])[minNumContributing] > .99
      #allZero = all(scaledWeights == 0)
    }
  }
  
  list(x = scaledWeights, ix = rawWeights$ix)
}

# Initialize the kernel at some small value based on the typical distances in the input distance matrix
minDistKernel = distMat[upper.tri(distMat)] %>% 
  unlist() %>% 
  sort() %>% #sort all observed distances
  .[. > 0] %>% 
  quantile(.001) # pick the .1th quantile. The only variants that will use this kernel will be in very densely populated regions of predictor space

varInfo %<>% 
  mutate(weightList = map(1:nrow(.), ~findWeights(.x, distMat, minDistKernel)))

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

## TODO: keep adapting over ulirschNegBinPrior
## 1. fit neg binomials
## 2. fit WEIGHTED gamma hyperprior to RNA counts
## 3. Fit marginal gamma hyperprior on DNA counts
## 4. adapt model code
## 5. run sampling
## 6. money
