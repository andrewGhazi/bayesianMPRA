#Let's see how much the informative prior helps using the output from allVarAnalysis.R

#ugh why do I do this crappy programming
varRes = varInference %>% sapply(function(x){x$hdi95}) %>% t %>% as.data.frame %>% as.tbl
names(varRes) = c('HDIlower', 'HDIupper')

varRes$construct = varInference %>% sapply(function(x){x$construct})
varRes$conclusion = varInference %>% sapply(function(x){x$bayesianConclusion})
varRes$wilcox.pval = varInference %>% sapply(function(x){x$wilcox.pval})
varRes$DeepSeaDnaase = varInference %>% sapply(function(x){x$constructData$DeepSeaDnaase})

#Restrict to those that were functional and in the same direction as DeepSea predicts. Let's just look through these manually
varRes %>% filter(conclusion == 'FUNCTIONAL') %>% filter(abs(DeepSeaDnaase) > .5, sign(HDIlower) == sign(DeepSeaDnaase))

