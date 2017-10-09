for (i in 1:7803) {
  fitGammaHyperPriors(i)
}

tmpfun = safely(fitGammaHyperPriors)
tmp = mclapply(1:7803, tmpfun, mc.cores = 20)

tmp = others %>% 
  mutate(weight = constr$weights[[1]]) %>% 
  dplyr::select(construct, negBinParams, weight) %>% 
  unnest %>% 
  na.omit %>% 
  filter(weight > 1e-4*sort(constr$weights[[1]], decreasing = TRUE)[30])

for (i in unique(tmp$type)) {
  for (j in unique(tmp$block)) {
    
    typeBlockDat = tmp %>%
      filter(type == i, block == j)
    
    # weightedLikFun = function(paramVec){-sum(typeBlockDat$weight * dgamma(typeBlockDat$muEst,
    #                                                   shape = paramVec[1],
    #                                                   rate = paramVec[2],
    #                                                   log = TRUE))}
    # 
    # plotDF = expand.grid(seq(.1, 10, length.out = 40),
    #                      seq(.01, .3, length.out = 40)) %>%
    #   as.tibble() %>%
    #   mutate(paramLik = map2_dbl(Var1, Var2, ~likFun(c(.x, .y))))
    # 
    # plotDF %>% ggplot(aes(Var1, Var2,)) + geom_raster(aes(fill = -paramLik))
    #   
    # typeBlockDat %$%
    #   fitdistMod(muEst, 
    #              'gamma', 
    #              start = initialMuGuess, weights = weight,
    #              control = list(ndeps = c(1e-4, 1e-5)),
    #              lower = 1e-10)
    
    typeBlockDat %$%
      fitdistMod(sizeEst[sizeEst < 1e4], 
                 'gamma', 
                 weights = weight[sizeEst < 1e4],
                 start = initialSizeGuess,
                 lower = 1e-10)
  }
}
