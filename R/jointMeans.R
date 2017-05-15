#Let's try to put together an example figure of Bayesian inference with the
#conditional density prior. So we have to actually do it. Basing this code off
#of MPRAMCMC.RMd in the MPRA rstudio project

library(dplyr)
library(mcmc) #Let's try to actually use this package this time instead of reinventing the wheel
library(parallel)
library(ggplot2)
library(coda)
filter = dplyr::filter
par(xpd = FALSE)

load("~/Qual/outputs/bandwidthOptimization9_27_16.RData")
load("~/MPRA/data/varInfo.RData")
load("~/MPRA/data/gatheredMPRA.RData")

h = which.max(resmat) %>% arrayInd(dim(resmat))
DSD.dens = density(varInfo$DeepSeaDnaase)

HPDinterval <- function(obj, prob = 0.95, ...) UseMethod("HPDinterval")

HPDinterval.mcmc <- function(obj, prob = 0.95, ...)
{
  #Adapted from the HPDinterval() function in the coda package
  obj <- obj %>% as.matrix()
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals)) stop("obj must have nsamp > 1")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                2, which.min)
  ans <- cbind(vals[cbind(inds, 1:npar)],
               vals[cbind(inds + gap, 1:npar)])
  dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
  attr(ans, "Probability") <- gap/nsamp
  ans
}

postfun <- function(ini, priorfun, y, type){
  #log-Posterior function
  #refMu - proposed reference mean
  #mutMu - proposed mutant mean
  #sigma - proposed variance
  #priorfun - after doing the conditional density estimation, use approxfun on the density to get this
  #y - vector of outcomes (quantile normalized activity levels)
  #type - factor with levels c('Ref', 'Mut'), corresponding to the type values of y
  #sig - standard deviation of activity after acounting for type
  
  refMu = ini[1]
  mutMu = ini[2]
  sigma = ini[3]
  
  ref.llik = sum(type == 'Ref') * log(sqrt(1/(2*pi*sigma^2))) + sum(-(y[type == 'Ref'] - refMu)^2/(2*sigma^2))
  mut.llik = sum(type == 'Mut') * log(sqrt(1/(2*pi*sigma^2))) + sum(-(y[type == 'Mut'] - mutMu)^2/(2*sigma^2))
  # llik = length(y)*log(sqrt(1/(2*pi*sig^2))) + sum(-(y - shiftvec - center)^2/(2*sig^2)) 
  
  TS = mutMu - refMu
  logPrior = log(priorfun(TS))
  
  return(ref.llik + mut.llik + logPrior) 
}

#Evaluate the likelihood
likefun <- function(ini, y, type, ordering){
  
  refMu = ini[1]
  mutMu = ini[2]
  sigma = ini[3]

  ref.llik = sum(type == 'Ref') * log(sqrt(1/(2*pi*sigma^2))) + sum(-(y[type == 'Ref'] - refMu)^2/(2*sigma^2))
  mut.llik = sum(type == 'Mut') * log(sqrt(1/(2*pi*sigma^2))) + sum(-(y[type == 'Mut'] - mutMu)^2/(2*sigma^2))
  
  return(ref.llik + mut.llik)
}


analyze <- function(constr, dirstr = 'varOutputs/plots/', tryMarg = FALSE, verbose = TRUE){ #, run = NULL, lik.run = NULL
  #Perform MCMC to estimate posterior ref/mut means based on empirical Bayesian prior. TS posterior from run$batch[,1] - run$batch[,2]
  #Also output wilcox p-value for comparison
  if(verbose){print('Setting up')}
  tmpVarInfo = varInfo %>% filter(construct != constr)
  testdat = varInfo %>% filter(construct == constr)
  mainparts = strsplit(constr, ' ') %>% unlist
  mainstr = paste0('chr', mainparts[1], ' position ', mainparts[2], ' ', testdat$ref, ' --> ', testdat$alt)
  
  #Let's see how the frequentist test does
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
  pfun = approxfun(prior.dens$x, prior.dens$y) #prior function
  
  qndat = MPRA.qnactivity %>% filter(construct == constr)
  
  nburn = 5000
  nrun = 10000
  thin = 10
  scaleval = .15*c(1,1,.5)
  lik.scale = .15*c(1,1,.5)
  
  if(verbose){print('Starting likelihood evaluation...')}
  lik.strt = Sys.time()
  lik.burn = metrop(likefun,
                    initial = c(-0.1321954, -0.1366325, 1.6), #Initialize at the global means
                    nbatch = nburn,
                    scale = lik.scale,
                    nspac = thin,
                    y = qndat$qnact,
                    type = qndat$type)
  # print(paste0('Likelihood burn done in ', format(Sys.time() - lik.strt)))
  
  lik.burntime = Sys.time() - lik.strt
  lik.runstrt = Sys.time()
  lik.run = metrop(lik.burn,
                   nbatch = nrun,
                   scale = lik.scale,
                   nspac = thin,
                   y = qndat$qnact,
                   type = qndat$type)
  
  lik.runTime = Sys.time() - lik.runstrt
  lik.tot = Sys.time() - lik.strt
  # print(paste0('Likelihood Run done in ', format(lik.runTime)))
  # print(paste0('Likelihood Total time ', format(lik.tot)))
  if(verbose){print('Done.')
    print('Starting posterior MCMC...')}
  #Evaluate the posterior
  strt = Sys.time()
  burn = metrop(postfun,
                initial = c(-0.1321954, -0.1366325, 1.6), #Initialize at the global means
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
  
  if(tryMarg){
    print('Starting posterior MCMC with marginal prior for comparison purposes...')
    marg.prior.dens = density(tmpVarInfo$VariantShift, from = -6, to = 6) #This is the optimized prior for the test construct!!!
    marg.pfun = approxfun(marg.prior.dens$x, marg.prior.dens$y) #marginal prior function
    margburn = metrop(postfun,
                      initial = c(-0.1321954, -0.1366325, 1.6), #Initialize at the global means
                      nbatch = nburn,
                      scale = scaleval,
                      nspac = thin,
                      y = qndat$qnact,
                      type = qndat$type,
                      priorfun = marg.pfun)
    margrun = metrop(margburn,
                     nbatch = nrun,
                     scale = scaleval,
                     nspac = thin,
                     y = qndat$qnact,
                     type = qndat$type,
                     priorfun = marg.pfun)
    margTSpost = margrun$batch[,2] - margrun$batch[,1]
    margHDI = HPDinterval(mcmc(margTSpost))
    
    if(margHDI[1] < 0 && margHDI[2] > 0){
      margConclusion = 'NONFUNCTIONAL'
    } else{
      margConclusion = 'FUNCTIONAL'
    }
    
    print('Done.')
  }
  
  tot = Sys.time() - strt
  # print(paste0('Posterior Run done in ', format(runTime)))
  # print(paste0('Posterior total time ', format(tot)))
  if(verbose){
    print('Done.')
    print('Starting output evaluation...')
  }
  
  
  tsPost = run$batch[,2] - run$batch[,1]
  tsLike = lik.run$batch[,2] - lik.run$batch[,1]
  
  hdi = HPDinterval(mcmc(tsPost))
  ymx = max(prior.dens$y, density(tsLike)$y, density(tsPost)$y)
  
  png(filename = paste0('outputs/', dirstr, (constr %>% gsub(' ', '_', .) %>% gsub('/', '-', .)), '.png'),
      width = 1280,
      height = 960)
  plot(density(tsPost), 
       col = 'forestgreen',
       lwd = 2,
       xlim = c(-2,2),
       ylim = c(0,ymx), main = mainstr,
       xlab = 'Transcriptional Shift',
       zero.line = FALSE,
       lend = 'butt')
  lines(density(tsLike),
        col = 'dodgerblue2',
        lwd = 2,
        lend = 'butt')
  lines(prior.dens$x,
        prior.dens$y,
        lwd = 2,
        col = 'firebrick2',
        lend = 'butt')
  lines(density(abs(run$batch[,3]))$x,
        .1*density(abs(run$batch[,3]))$y,
        col = 'chocolate1',
        lwd = 2,
        lend = 'butt')
  lines(hdi[1,],
        rep(-.02235 * ymx,2),
        col = 'darkorchid2',
        lwd = 4,
        lend = 'butt')
  legend('topright',
         legend =  c('Prior', 'Likelihood', 'Posterior', '95% HDI', 'Activity Std Dev'),
         fill = c('firebrick2', 'dodgerblue2', 'forestgreen', 'darkorchid2', 'chocolate1'),
         cex = 1.25)
  dev.off()
  
  if(hdi[1] < 0 && hdi[2] > 0){
    
    bayes.call = 'NONFUNCTIONAL'
  } else{
    
    bayes.call = 'FUNCTIONAL'
  }
  
  
  if(tryMarg){
    res = list(constr, testdat, run, prior.dens, lik.run, hdi, bayes.call, margrun, margHDI, margConclusion, pwilcox)
    names(res) = c('construct', 'constructData','posterior', 'prior', 'likelihood', 'hdi95', 'bayesianConclusion', 'posteriorWithMarginalPrior', 'margHDI', 'margConclusion', 'wilcox.pval')
  } else{
    res = list(constr, testdat, run, prior.dens, lik.run, hdi, bayes.call, pwilcox)
    names(res) = c('construct', 'constructData','posterior', 'prior', 'likelihood', 'hdi95', 'bayesianConclusion', 'wilcox.pval')
  }
  
  if(verbose){
    print('All Done!')
  }
  return(res)
}

# this function deprecated. Use tryMarg = TRUE in analyze()
# marginalPriorAnalyze <- function(constr, dirstr = 'varOutputs/plots/', verbose = TRUE){ #, run = NULL, lik.run = NULL
#   #Perform MCMC to estimate posterior ref/mut means based on empirical Bayesian prior. TS posterior from run$batch[,1] - run$batch[,2]
#   #Also output wilcox p-value for comparison
#   if(verbose){print('Setting up')}
#   tmpVarInfo = varInfo %>% filter(construct != constr)
#   testdat = varInfo %>% filter(construct == constr)
#   mainparts = strsplit(constr, ' ') %>% unlist
#   mainstr = paste0('chr', mainparts[1], ' position ', mainparts[2], ' ', testdat$ref, ' --> ', testdat$alt)
#   
#   #Let's see how the frequentist test does
#   observations = MPRA.qnactivity %>% filter(construct == constr)
#   #t.test(qnact~type, data = observations) # p = 0.000777 -- It doesn't get below the significance threshold used by Ulirsch (9e-5)
#   pwilcox = wilcox.test(observations %>% filter(type == 'Ref') %>% .$qnact,
#                         observations %>% filter(type == 'Mut') %>% .$qnact)$p.value
#   #stripchart(qnact ~ type, data = observations, vertical = TRUE, pch =16, method = 'jitter', main = mainstr)
#   
#   ####### Doing the MCMC
#   
#   
#   inv.sqrt.dens = approx(DSD.dens$x, 1/((DSD.dens$y)^(1/3)), xout = testdat$DeepSeaDnaase)$y
#   w = dnorm(tmpVarInfo$DeepSeaDnaase, mean = testdat$DeepSeaDnaase, sd = h2.grid[h[2]]*inv.sqrt.dens) #incorporating the condition location = adaptive
#   prior.dens = density(tmpVarInfo$VariantShift, from = -6, to = 6) #This is the optimized prior for the test construct!!!
#   # , weights = w/sum(w), kernel = 'gaussian', bw = h1.grid[h[1]], cut this out from analyze(). Now the prior is just the marginal distribution
#   # plot(prior.dens,
#   #      main = paste0('Prior on ', mainstr))
#   pfun = approxfun(prior.dens$x, prior.dens$y)
#   
#   qndat = MPRA.qnactivity %>% filter(construct == constr)
#   
#   nburn = 5000
#   nrun = 10000
#   thin = 10
#   scaleval = .15*c(1,1,.5)
#   lik.scale = .2
#   
#   if(verbose){print('Starting likelihood evaluation...')}
#   lik.strt = Sys.time()
#   lik.burn = metrop(likefun,
#                     initial = c(-0.1321954, -0.1366325, 1.6),
#                     nbatch = nburn,
#                     scale = lik.scale,
#                     nspac = thin,
#                     y = qndat$qnact,
#                     type = qndat$type)
#   # print(paste0('Likelihood burn done in ', format(Sys.time() - lik.strt)))
#   
#   lik.burntime = Sys.time() - lik.strt
#   lik.runstrt = Sys.time()
#   lik.run = metrop(lik.burn,
#                    nbatch = nrun,
#                    scale = lik.scale,
#                    nspac = thin,
#                    y = qndat$qnact,
#                    type = qndat$type)
#   
#   lik.runTime = Sys.time() - lik.runstrt
#   lik.tot = Sys.time() - lik.strt
#   # print(paste0('Likelihood Run done in ', format(lik.runTime)))
#   # print(paste0('Likelihood Total time ', format(lik.tot)))
#   if(verbose){print('Done.')
#     print('Starting posterior MCMC...')}
#   #Evaluate the posterior
#   strt = Sys.time()
#   burn = metrop(postfun,
#                 initial = c(-0.1321954, -0.1366325, 1.6),
#                 nbatch = nburn,
#                 scale = scaleval,
#                 nspac = thin,
#                 y = qndat$qnact,
#                 type = qndat$type,
#                 priorfun = pfun)
#   # print(paste0('Posterior Burn done in ', format(Sys.time() - strt)))
#   
#   burntime = Sys.time() - strt
#   runstrt = Sys.time()
#   run = metrop(burn,
#                nbatch = nrun,
#                scale = scaleval,
#                nspac = thin,
#                y = qndat$qnact,
#                type = qndat$type,
#                priorfun = pfun)
#   runTime = Sys.time() - runstrt
#   
#   tot = Sys.time() - strt
#   # print(paste0('Posterior Run done in ', format(runTime)))
#   # print(paste0('Posterior total time ', format(tot)))
#   if(verbose){
#     print('Done.')
#     print('Starting output evaluation...')
#   }
#   
#   
#   tsPost = run$batch[,2] - run$batch[,1]
#   tsLike = lik.run$batch[,2] - lik.run$batch[,1]
#   
#   hdi = HPDinterval(mcmc(tsPost))
#   ymx = max(prior.dens$y, density(tsLike)$y, density(tsPost)$y)
#   
#   png(filename = paste0('outputs/', dirstr, 'margPrior', (constr %>% gsub(' ', '_', .) %>% gsub('/', '-', .)), '.png'),
#       width = 1280,
#       height = 960)
#   plot(density(tsPost), 
#        col = 'forestgreen',
#        lwd = 2,
#        xlim = c(-2,2),
#        ylim = c(0,ymx), main = mainstr,
#        xlab = 'Transcriptional Shift',
#        zero.line = FALSE,
#        lend = 'butt')
#   lines(density(tsLike),
#         col = 'dodgerblue2',
#         lwd = 2,
#         lend = 'butt')
#   lines(prior.dens$x,
#         prior.dens$y,
#         lwd = 2,
#         col = 'firebrick2',
#         lend = 'butt')
#   lines(density(abs(run$batch[,3]))$x,
#         .1*density(abs(run$batch[,3]))$y,
#         col = 'chocolate1',
#         lwd = 2,
#         lend = 'butt')
#   lines(hdi[1,],
#         rep(-.02235 * ymx,2),
#         col = 'darkorchid2',
#         lwd = 4,
#         lend = 'butt')
#   legend('topright',
#          legend =  c('Prior', 'Likelihood', 'Posterior', '95% HDI', 'Activity Std Dev'),
#          fill = c('firebrick2', 'dodgerblue2', 'forestgreen', 'darkorchid2', 'chocolate1'),
#          cex = 1.25)
#   dev.off()
#   
#   if(hdi[1] < 0 && hdi[2] > 0){
#     
#     bayes.call = 'NONFUNCTIONAL'
#   } else{
#     
#     bayes.call = 'FUNCTIONAL'
#   }
#   
#   
#   
#   res = list(constr, testdat, run, prior.dens, lik.run, hdi, bayes.call, pwilcox)
#   names(res) = c('construct', 'constructData','posterior', 'prior', 'likelihood', 'hdi95', 'bayesianConclusion', 'wilcox.pval')
#   if(verbose){
#     print('All Done!')
#   }
#   return(res)
# }