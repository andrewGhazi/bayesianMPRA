## Let's try to do the predictor variable selection analysis thing
## Let's see what subset of predictors performs the best by the conditional/marginal density ratio metric

library(glmnet)
library(tidyverse)

load("~/bayesianMPRA/analysis_data/varInfoWithHistoneMarkAnnotations.RData")

model_data = varInfo %>% 
  select(transcription_shift, eigen:DeepSeaDnaase, matches('Broad|Syd')) %>% 
  na.omit

predictors = names(model_data)[-1]

lasso_model = cv.glmnet(alpha = 1,
                        x = model_data %>% select(-transcription_shift) %>% as.matrix,
                        y = model_data %>% select(transcription_shift) %>% as.matrix,
                        nfolds = 10)

plot(lasso_model,
     label = TRUE)
