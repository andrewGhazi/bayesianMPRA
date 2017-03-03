permtest = vector('numeric', 100000)
refdat = observations %>% filter(type == 'Ref')
mutdat = observations %>% filter(type == 'Mut')
for(i in 1:100000){
  permtest[i] = mean(base::sample(mutdat$qnact, 140, replace = TRUE)) - mean(base::sample(refdat$qnact, 140, replace = TRUE))
}
hist(permtest, breaks = 100)

postfun <- function(TSsig, priorfun, y, type){
  #log-Posterior function
  #TS - proposed transcriptional shift
  #priorfun - after doing the conditional density estimation, use approxfun on the density to get this
  #y - vector of outcomes (quantile normalized activity levels)
  #type - factor with levels c('Ref', 'Mut'), corresponding to the type values of y
  #sig - standard deviation of activity after acounting for type
  #mult - multiplier to make exp(dEigen) * exp(dGkmer) be a probability distribution
  
  
  #sig = sd(y[type == 'Ref'])
  TS = TSsig[1]
  sig = TSsig[2]
  
  shiftvec = rep(TS, length(y))
  shiftvec[type == 'Ref'] = rep(0, sum(type == 'Ref'))
  center = mean(y[type == 'Ref']) #TODO: do this in an unbiased way
  llik = length(y)*log(sqrt(1/(2*pi*sig^2))) + sum(-(y - shiftvec - center)^2/(2*sig^2)) 
  
  
  logPrior = log(priorfun(TS))
  
  return(llik + logPrior) 
}

#Evaluate the likelihood
likefun <- function(TS, y, type, sig = 2.5, ordering){
  sig = sd(y[type == 'Ref'])
  
  shiftvec = rep(TS, length(y))
  shiftvec[type == 'Ref'] = rep(0, sum(type == 'Ref'))
  center = mean(y[type == 'Ref']) #TODO: do this in an unbiased way
  
  llik = length(y)*log(sqrt(1/(2*pi*sig^2))) + sum(-(y - shiftvec - center)^2/(2*sig^2))
  
  return(llik)
}


analyze <- function(constr, dirstr = 'varOutputs/plots/'){ #, run = NULL, lik.run = NULL
  #Perform MCMC to estimate posterior transcriptional shift based on empirical Bayesian prior
  #Also output wilcox p-value for comparison
  
  tmpVarInfo = varInfo %>% filter(construct != constr)
  testdat = varInfo %>% filter(construct == constr)
  mainparts = strsplit(constr, ' ') %>% unlist
  mainstr = paste0('chr', mainparts[1], ' position ', mainparts[2], ' ', testdat$ref, ' --> ', testdat$alt)
  
  #Let's see how the t.test does
  observations = MPRA.qnactivity %>% filter(construct == constr)
  #t.test(qnact~type, data = observations) # p = 0.000777 -- It doesn't get below the significance threshold used by Ulirsch (9e-5)
  pwilcox = wilcox.test(observations %>% filter(type == 'Ref') %>% .$qnact,
                        observations %>% filter(type == 'Mut') %>% .$qnact)$p.value
  #stripchart(qnact ~ type, data = observations, vertical = TRUE, pch =16, method = 'jitter', main = mainstr)
  
  ####### Doing the MCMC
  
  
  inv.sqrt.dens = approx(DSD.dens$x, 1/((DSD.dens$y)^(1/3)), xout = testdat$DeepSeaDnaase)$y
  w = dnorm(tmpVarInfo$DeepSeaDnaase, mean = testdat$DeepSeaDnaase, sd = h2.grid[h[2]]*inv.sqrt.dens) #incorporating the condition location = adaptive
  prior.dens = density(tmpVarInfo$VariantShift, weights = w/sum(w), kernel = 'gaussian', bw = h1.grid[h[1]], from = -6, to = 6) #This is the optimized prior for the test construct!!!
  # plot(prior.dens,
  #      main = paste0('Prior on ', mainstr))
  pfun = approxfun(prior.dens$x, prior.dens$y)
  
  
  
  qndat = MPRA.qnactivity %>% filter(construct == constr)
  
  nburn = 5000
  nrun = 10000
  thin = 10
  scaleval = .15*c(1,.5)
  
  lik.strt = Sys.time()
  lik.burn = metrop(likefun,
                    initial = 0,
                    nbatch = nburn,
                    scale = .8,
                    nspac = thin,
                    y = qndat$qnact,
                    type = qndat$type)
  # print(paste0('Likelihood burn done in ', format(Sys.time() - lik.strt)))
  
  lik.burntime = Sys.time() - lik.strt
  lik.runstrt = Sys.time()
  lik.run = metrop(lik.burn,
                   nbatch = nrun,
                   scale = .8,
                   nspac = thin,
                   y = qndat$qnact,
                   type = qndat$type)
  
  lik.runTime = Sys.time() - lik.runstrt
  lik.tot = Sys.time() - lik.strt
  # print(paste0('Likelihood Run done in ', format(lik.runTime)))
  # print(paste0('Likelihood Total time ', format(lik.tot)))
  
  #Evaluate the posterior
  strt = Sys.time()
  burn = metrop(postfun,
                initial = c(0,1.6),
                nbatch = nburn,
                scale = scaleval,
                nspac = thin,
                y = qndat$qnact,
                type = qndat$type,
                priorfun = pfun)
  # print(paste0('Posterior Burn done in ', format(Sys.time() - strt)))
  
  burntime = Sys.time() - strt
  runstrt = Sys.time()
  run = metrop(burn,
               nbatch = nrun,
               scale = scaleval,
               nspac = thin,
               y = qndat$qnact,
               type = qndat$type,
               priorfun = pfun)
  runTime = Sys.time() - runstrt
  tot = Sys.time() - strt
  # print(paste0('Posterior Run done in ', format(runTime)))
  # print(paste0('Posterior total time ', format(tot)))
  
  
  
  hdi = HPDinterval(run)
  
  ymx = max(prior.dens$y, density(lik.run$batch)$y, density(run$batch)$y)
  
  png(filename = paste0('outputs/', dirstr, (constr %>% gsub(' ', '_', .) %>% gsub('/', '-', .)), '.png'),
      width = 960,
      height = 960)
  plot(density(run$batch[,1]), 
       col = 'forestgreen',
       lwd = 2,
       xlim = c(-2,2),
       ylim = c(0,ymx), main = mainstr,
       xlab = 'Transcriptional Shift',
       zero.line = FALSE)
  lines(density(lik.run$batch),
        col = 'dodgerblue2',
        lwd = 2)
  lines(prior.dens$x,
        prior.dens$y,
        lwd = 2,
        col = 'firebrick2')
  lines(density(abs(run$batch[,2]))$x,
        .1*density(abs(run$batch[,2]))$y,
        col = 'chocolate1',
        lwd = 2)
  lines(hdi[1,],
        rep(-.05,2),
        col = 'darkorchid2',
        lwd = 4,
        lend = 1)
  legend('topright',
         legend =  c('Prior', 'Likelihood', 'Posterior', '95% HDI'),
         fill = c('firebrick2', 'dodgerblue2', 'forestgreen', 'darkorchid2'),
         cex = .85)
  dev.off()
  
  if(hdi[1] < 0 && hdi[2] > 0){
    
    bayes.call = 'NONFUNCTIONAL'
  } else{
    
    bayes.call = 'FUNCTIONAL'
  }
  
  
  
  res = list(constr, testdat, run, prior.dens, lik.run, hdi, bayes.call, pwilcox)
  names(res) = c('Construct', 'ConstructData','Posterior', 'Prior', 'Likelihood', 'HDI95', 'BayesianConclusion', 'Wilcox.pval')
  return(res)
  
}