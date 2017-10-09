library(dplyr)
library(magrittr)
library(parallel)
#source('BayesianAnalysis.R')
source('scripts/jointMeans.R')

load("~/MPRA/data/varInfo.RData")
load("~/MPRA/data/gatheredMPRA.RData")

varInference = mclapply(varInfo$construct %>% as.character, function(x){
  return(analyze(x, tryMarg = TRUE, verbose = FALSE))
}, mc.cores = 20)
save(varInference, 
     file = 'outputs/varOutputs/distributions/UlirschVariantJointMeansWithMarginal_Analysis.RData')


sapply(varInference, function(x){x$bayesianConclusion}) %>% table
# .
# FUNCTIONAL NONFUNCTIONAL 
# 2419          4950 


priorSource = data_frame(construct = sapply(varInference, function(x){x$Construct}),
                         ConditionalConclusion = sapply(varInference, function(x){x$bayesianConclusion}),
                         MarginalConclusion = sapply(varInference, function(x){x$margConclusion}))

table(priorSource[,2:3])

res = data_frame(construct = sapply(varInference, function(x){x$construct}),
                 ConditionalConclusion = sapply(varInference, function(x){x$bayesianConclusion}),
                 MarginalConclusion = sapply(varInference, function(x){x$margConclusion}),
                 wilcox.pval = sapply(varInference, function(x){x$wilcox.pval}))

res %<>% mutate(qValue = p.adjust(wilcox.pval, method = 'fdr'),
                freqConclusion95 = (qValue < .05) %>% ifelse(., 'FUNCTIONAL', 'NONFUNCTIONAL'),
                freqConclusion99 = (qValue < .01) %>% ifelse(., 'FUNCTIONAL', 'NONFUNCTIONAL'))#FINE HADLEY I ADMIT THIS DPLYR STUFF IS NICE WHEN YOU GET THE HANG OF IT

table(res %>% select(ConditionalConclusion, freqConclusion95)) #Comparison of frequentist and empirical Bayesian calls

#Let's look at the 99% HDI's to compare to Ulirsch's analysis
hdi99 = sapply(varInference, function(x){HPDinterval(mcmc(x$posterior$batch[,2] - x$posterior$batch[,1]), prob = .99)}) %>% t

res$hdi99lower = hdi99[,1]
res$hdi99upper = hdi99[,2]
res %<>% mutate(hdi99Conclusion = ifelse(hdi99lower < 0 & hdi99upper > 0, 'NONFUNCTIONAL', 'FUNCTIONAL'))

table(res %>% select(hdi99Conclusion, freqConclusion99))

#Let's compute posterior error probabilities as per http://varianceexplained.org/r/bayesian_fdr_baseball/
# a PEP is in this case the probability of committing a type S error (calling a variant as functional when the true sign of the variant is on the other side of 0)
res$PEP = mclapply(varInference, #Posterior Error Probability
                   function(x){
                     if(mean(x$posterior$batch[,2] - x$posterior$batch[,1]) < 0){ #if it's an activity decreasing construct
                       return(1 - ecdf(x$posterior$batch[,2] - x$posterior$batch[,1])(0)) #evaluate the empirical distribution function for TS > 0
                     }else{ #if it's an activity increasing construct
                       return(ecdf(x$posterior$batch[,2] - x$posterior$batch[,1])(0)) # evaluate the ecdf for TS < 0
                     }
                   },
                   mc.cores = 20) %>% unlist

#res %<>% arrange(PEP) %>% mutate(bayesianFDR = cummean(PEP))
#res %>% filter(bayesianFDR < .01) #1175 constructs, so this is very generous. Makes sense given that the criteria for inclusion is just that the sign is correct rather than being far from 0.

margHDI99 = sapply(varInference,
                     function(x){
                       return(HPDinterval(mcmc(x$posteriorWithMarginalPrior$batch[,2] - x$posteriorWithMarginalPrior$batch[,1]), prob = .99))
                     }) %>% t
res$margHDI99lower = margHDI99[,1]
res$margHDI99upper = margHDI99[,2]
res %<>% mutate(marg99Conclusion = ifelse(margHDI99lower < 0 & margHDI99upper > 0, 'NONFUNCTIONAL', 'FUNCTIONAL'))

res$DSD = sapply(varInference,#DeepSeaDNase I hypersensitivity (predicted)
                 function(x){x$constructData$DeepSeaDnaase[1]})
res$TS = sapply(varInference, #Transcriptional shift
                 function(x){x$constructData$VariantShift[1]})

#This shows the (76) constructs that get called by a 99%HDI with the conditional empirical prior that don't get called under any other circumstances and show their stats
# Hopefully adding more predictors to the conditional prior will produce more of these.
par(mar = c(4.5,5,4,1) + .1)
res %>% filter(marg99Conclusion == 'NONFUNCTIONAL' & hdi99Conclusion == 'FUNCTIONAL' & freqConclusion99 == 'NONFUNCTIONAL') %>% # & abs(DSD) > .5
  select(construct, hdi99lower, DSD, TS) %T>% 
  plot(TS ~ DSD, data = .,
       xlab = 'DeepSea predicted DNase I hypersensitivity',
       ylab = 'Observed mean activity difference\n(Transcriptional Shift)') %>% 
  summary(lm(TS ~ DSD, data = .))
