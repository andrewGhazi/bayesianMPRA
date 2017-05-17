library(tidyverse)
library(magrittr)
library(stringr)

load('data/varInfoWithHistoneMarkAnnotations.RData')
names(varInfo)[12] = 'transcriptionalShift'

preds = varInfo %>% select(contains('Broad'), contains('Sydh')) %>%
  map_df(~scale(.x)[,1]) #effing scale only outputs matrices

generateDistMat = function(predictors) {
  #predictors is a n x d data frame of predictors
  
  predictors %>% 
    as.data.frame() %>% 
    as.matrix %>% 
    dist %>% 
    as.matrix
  
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

varInfo %<>% 
  mutate(kNN = map(1:nrow(.), ~findKNN(k, .x, distMat)),
         numNeigh = map_int(kNN, length))

computeLikelihood

