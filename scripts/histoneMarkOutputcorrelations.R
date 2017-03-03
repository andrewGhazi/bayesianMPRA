library(purrr)
library(dplyr)
library(magrittr)
load("~/Qual/data/varInfoWithHistoneMarkAnnotations.RData")

varInfo %>% 
  select(11:27) %>% 
  select(-VariantShift) %>% 
  gather(markType, value, -DeepSeaDnaase) %>% 
  split(.$markType) %>% 
  map(~ lm(DeepSeaDnaase ~ value, data = .)) %>% 
  map(summary)

varInfo %>% 
  select(12:27) %>% 
  gather(markType, value, -VariantShift) %>% 
  split(.$markType) %>% 
  map(~ lm(VariantShift ~ value, data = .)) %>% 
  map(summary)

varInfo %>% 
  select(VariantShift, BroadK562H3k36me3) %>% 
  ggplot(aes(BroadK562H3k36me3, VariantShift)) +
  geom_point()

#So in isolation none of these are amazing predictors for transcriptional shift
# Let's try LASSO
library(glmnet)
dat = varInfo %>% 
  select(12:27) %>% as.matrix
fit = glmnet(dat[,-1], dat[,1],
             alpha = 1)

cv.fit = cv.glmnet(dat[,-1], dat[,1],
                   alpha = 1)
