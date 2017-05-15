library(dplyr)
library(ggplot2)
library(magrittr)
load("~/MPRA/data/gatheredMPRA.RData")


#tmp = analyze('11 8902218 1/2', dirstr = 'tmp/')
constr = '17 37520449 1/2'
#constr = '6 32609097 1/3'

MPRA.qnactivity %>% #data plot
  filter(construct == constr) %>% 
  ggplot(data = ., mapping = aes(type, qnact)) + 
  geom_violin() + 
  geom_jitter(width = .1, height = 0)

tmpVarInfo = varInfo %>% filter(construct != constr)
testdat = varInfo %>% filter(construct == constr)
tmp = varInference[[which(sapply(varInference, function(x){x$construct}) == constr)]]
run = tmp$posterior
lik.run = tmp$likelihood
prior.dens = tmp$prior

tsPost = run$batch[,2] - run$batch[,1]
tsLike = lik.run$batch[,2] - lik.run$batch[,1]

hdi = HPDinterval(mcmc(tsPost), prob = .99)
ymx = max(prior.dens$y, density(tsLike)$y, density(tsPost)$y)

mainparts = strsplit(constr, ' ') %>% unlist
mainstr = paste0('chr', mainparts[1], ' position ', mainparts[2], ' ', testdat$ref, '     ', testdat$alt)


png(filename = 'outputs/figures/example.png',
   width = 1080,
   height = 810,
   res = 150)
par(mgp = c(2.5, 1, 0))
plot(density(tsPost), 
     col = 'forestgreen',
     lwd = 2,
     xlim = c(-2,2),
     ylim = c(0,ymx),
     main = expression(chr17 ~ position ~ 37520449 ~ T %->% C),
     xlab = 'Transcriptional Shift',
     zero.line = FALSE,
     lend = 'butt',
     cex.lab = 1.5,
     cex.main = 1.5)
lines(density(tsLike),
      col = 'dodgerblue2',
      lwd = 2,
      lend = 'butt')
lines(prior.dens$x,
      prior.dens$y,
      lwd = 2,
      col = 'firebrick2',
      lend = 'butt')
# lines(density(abs(run$batch[,3]))$x,
#       .1*density(abs(run$batch[,3]))$y,
#       col = 'chocolate1',
#       lwd = 2,
#       lend = 'butt')
lines(hdi[1,],
      rep(-.02235 * ymx,2),
      col = 'darkorchid2',
      lwd = 4,
      lend = 'butt')
legend('topright',
       legend =  c('Prior', 'Likelihood', 'Posterior', '99% HDI'), #, 'Activity Std Dev'
       fill = c('firebrick2', 'dodgerblue2', 'forestgreen', 'darkorchid2'),#, 'chocolate1'
       cex = 1.25)
dev.off()

#MFN2 expression association
library(dplyr)
library(magrittr)
load('/mnt/labhome/simonlm/projects/PRAX/Papers/ResourcePaper/data/mRNAsAboveCutoff.RData') #PRAX expression data

mRNA.info %<>% as.tbl
mfn2TranscriptID = mRNA.info %>% filter(Symbol == 'MFN2') %>% .$TranscriptClusterID
mfn2Row = which((mRNA.cutoff.mains %>% row.names) == mfn2TranscriptID)

rs1474868genotype = system('grep \'rs1474868\' /mnt/bigData/simonlm/Bray/PRAX/data/genotype/Consolidated.snps.txt', intern = TRUE)[1] %>% 
  strsplit(split = '\t') %>% 
  unlist
genotype = rs1474868genotype[2:155]

expr = mRNA.cutoff.mains[mfn2Row,]

png(filename = 'outputs/figures/MFN2plot.png',
    width = 1080,
    height = 810,
    res = 200)
par(mar = c(4,5,3,1)+.1)
stripchart(expr ~ genotype,
           group.names = c('CC', 'CT', 'TT'),
           vertical = TRUE,
           method = 'jitter',
           pch = 19,
           ylab = 'MFN2 Expression',
           xlab = 'Genotype',
           main = 'rs1474868',
           cex.lab = 1.5,
           cex = .8,
           jitter = .15)
dev.off()


attach(varInfo)

p = ggplot(varInfo, aes(DeepSeaDnaase, VariantShift)) + 
  geom_point()

cvals = seq(range(DeepSeaDnaase)[1], range(DeepSeaDnaase)[2], length.out = 10) #values to condition on -1:8
DSD.dens = density(DeepSeaDnaase)
inv.sqrt.dens = approx(DSD.dens$x, 1/sqrt(DSD.dens$y), xout = cvals)$y

for(i in 1:length(cvals)){
  w = dnorm(DeepSeaDnaase, mean = cvals[i], sd = h2.grid[h[2]]*inv.sqrt.dens[i]) #incorporating the condition location = adaptive
  dens = density(VariantShift, weights = w/sum(w), kernel = 'gaussian', bw = h1.grid[h[1]])
  if(i == 1){
    cdens = tibble(x = dens$x, y = dens$y, posval = rep(cvals[i], length(dens$x)))
  } else{
    cdens = rbind(cdens, tibble(x = dens$x, y = dens$y, posval = rep(cvals[i], length(dens$x))))
  }
}

p = p + geom_path(data = cdens, aes(x = .75*sd(DeepSeaDnaase)*y + posval, y = x, group = posval), color = 'firebrick1')
# wex = tibble(x = seq(range(DeepSeaDnaase)[1], range(DeepSeaDnaase)[2], length.out = 200),
#              y = .1*dnorm(x, mean = mean(range(DeepSeaDnaase)), sd = diff(range(DeepSeaDnaase))/ 30)) # weight example
#p = p + geom_line(data = wex, aes(x,y-3))
p + ylab('Transcriptional Shift') + xlab('DeepSea DNaseI Hypersensitivity') +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))
ggsave('conditionalDensBigLabels.png', dpi = 200)
detach(varInfo)