# Now let's try to run the hierarchical model but estimate the parameters of the prior from the 30kNN

library(tidyverse)
library(stringr)
library(magrittr)

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

neighbors = varInfo %>% 
  filter(construct == '1 155271258 1/3') %>% 
  .$kNN %>% 
  unlist

varInfo[neighbors,]$countData %>% #estimate mean/variances by neighbor, estimate Gamma distributions off of those

