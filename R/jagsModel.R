# Could we just define a JAGS model that works on the raw counts using a negative binomial model rather than the activity based stuff
# This would also naturally account for variants with barcodes at count 0

# NegBin reparameterization adapted from here: http://doingbayesiandataanalysis.blogspot.com/2012/04/negative-binomial-reparameterization.html

library(arm)
library(rjags)
library(coda)

set.seed(1280)

dir="/mnt/labhome/andrew/MPRA/paper_data/"

exampleData = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                         delim = "\t",
                         col_names = T,
                         col_types = cols(chr = "c")) %>% 
  dplyr::select(matches('construct|CTRL|DNA|type')) %>% 
  filter(construct == '1 155271258 1/3') %>% 
  mutate(bcIndex = 1:nrow(.)) 

dataList = list(nRefBarcode = sum(exampleData$type == 'Ref'),
                nMutBarcode = sum(exampleData$type == 'Mut'),
                nDNAblocks = colnames(exampleData) %>% grepl('DNA', .) %>% sum,
                nRNAblocks = colnames(exampleData) %>% grepl('RNA', .) %>% sum,
                refDNAmat = exampleData %>% filter(type == 'Ref') %>% select(contains('DNA')) %>% as.data.frame %>% as.matrix,
                refRNAmat = exampleData %>% filter(type == 'Ref') %>% select(contains('RNA')) %>% as.data.frame %>% as.matrix,
                mutDNAmat = exampleData %>% filter(type == 'Mut') %>% select(contains('DNA')) %>% as.data.frame %>% as.matrix,
                mutRNAmat = exampleData %>% filter(type == 'Mut') %>% select(contains('RNA')) %>% as.data.frame %>% as.matrix)


#--------------
modelStringInitial = "
model {
  for ( i in 1:nRefBarcode ) {
    for ( j in 1:nDNAblocks ) {
      refDNAmat[i,j] ~ dnegbin(pRefDNA, rRefDNA)
    }
    for ( j in 1:nRNAblocks ) {
      refRNAmat[i,j] ~ dnegbin(pRefRNA, rRefRNA)
    }
  }

  for ( i in 1:nMutBarcode ) {
    for ( j in 1:nDNAblocks ) {
      mutDNAmat[i,j] ~ dnegbin(pMutDNA, rMutDNA)
    }
    for ( j in 1:nRNAblocks ) {
      mutRNAmat[i,j] ~ dnegbin(pMutRNA, rMutRNA)
    }
  }

  pRefDNA <- rRefDNA/(rRefDNA+mRefDNA)
  rRefDNA ~ dgamma(0.01, 0.01)
  mRefDNA ~ dgamma(0.01, 0.01)
  vRefDNA <- rRefDNA*(1-pRefDNA)/(pRefDNA*pRefDNA)

  pRefRNA <- rRefRNA/(rRefRNA+mRefRNA)
  rRefRNA ~ dgamma(0.01, 0.01)
  mRefRNA ~ dgamma(0.01, 0.01)
  vRefRNA <- rRefRNA*(1-pRefRNA)/(pRefRNA*pRefRNA)

  pMutDNA <- rMutDNA/(rMutDNA+mMutDNA)
  rMutDNA ~ dgamma(0.01, 0.01)
  mMutDNA ~ dgamma(0.01, 0.01)
  vMutDNA <- rMutDNA*(1-pMutDNA)/(pMutDNA*pMutDNA)

  pMutRNA <- rMutRNA/(rMutRNA+mMutRNA)
  rMutRNA ~ dgamma(0.01, 0.01)
  mMutRNA ~ dgamma(0.01, 0.01)
  vMutRNA <- rMutRNA*(1-pMutRNA)/(pMutRNA*pMutRNA)
}
"

#---------------

jagsModel = jags.model(file = textConnection(modelStringInitial),
                       data = dataList,
                       n.chains = 3, 
                       n.adapt = 500)
update(jagsModel,
       n.iter = 500)
codaSamples = coda.samples(jagsModel,
                           variable.names = c('mRefDNA', 'mRefRNA', 'vRefDNA', 'vRefRNA', 'mMutDNA', 'mMutRNA', 'vMutDNA', 'vMutRNA'),
                           n.iter = 3334)

res = codaSamples %>% 
  map_df(~as.matrix(.x) %>% as.data.frame %>% as.tbl) %>% 
  mutate(TS = log(mMutRNA / mMutDNA) - log(mRefRNA / mRefDNA))

res %>% 
  ggplot(aes(TS)) + geom_histogram(bins = 50)

# Hierarchical model ------------
# Hierarchical because model parameters come from a common distribution
# Probably doesn't make sense to have the DNA / RNA come from the same distribution but it's a start
modelStringHier = "
model {
  for ( i in 1:nRefBarcode ) {
    for ( j in 1:nDNAblocks ) {
      refDNAmat[i,j] ~ dnegbin(p[1], r[1])
    }
    for ( j in 1:nRNAblocks ) {
      refRNAmat[i,j] ~ dnegbin(p[2], r[2])
    }
  }

  for ( i in 1:nMutBarcode ) {
    for ( j in 1:nDNAblocks ) {
      mutDNAmat[i,j] ~ dnegbin(p[3], r[3])
    }
    for ( j in 1:nRNAblocks ) {
      mutRNAmat[i,j] ~ dnegbin(p[4], r[4])
    }
  }

  for ( i in 1:4 ) { #1 = refDNA, 2 = refRNA, 3 = mutDNA, 4 = mutRNA
    p[i] <- r[i]/(r[i] + m[i])
    r[i] ~ dgamma(0.01, 0.01)
    m[i] ~ dgamma(0.01, 0.01)
    v[i] <- r[i]*(1-p[i])/(p[i]*p[i])
  }
}
"

#---------------

jagsModelHier = jags.model(file = textConnection(modelStringHier),
                           data = dataList,
                           n.chains = 3, 
                           n.adapt = 500)
update(jagsModelHier,
       n.iter = 500)
codaSamplesHier = coda.samples(jagsModelHier,
                               variable.names = c('m', 'v'),
                               n.iter = 3334)

res = codaSamples %>% 
  map_df(~as.matrix(.x) %>% as.data.frame %>% as.tbl) %>% 
  mutate(TS = log(mMutRNA / mMutDNA) - log(mRefRNA / mRefDNA))

res %>% 
  ggplot(aes(TS)) + geom_histogram(bins = 50)