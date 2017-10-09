tmp = '6 109625879 2/3'
tmpVar = varFuns %>% filter(construct == tmp) 
likFun = tmpVar$varLikFun[[1]]
priorFun = tmpVar$priorFun[[1]]
postFun = tmpVar$postFun[[1]]

data_frame(x = seq(-10, 10, by = .1), 
           lik = map_dbl(x, ~likFun(c(.x, 1, 1.5))), 
           prior = map_dbl(x, ~priorFun(c(.x, 1, 1.5))), 
           post = map_dbl(x, ~postFun(c(.x, 1, 1.5)))) %>% 
  gather(dist, dens, -x) %>% 
  ggplot(aes(x, dens, color = dist)) + 
  geom_line() + 
  ylab('log density') + 
  xlab('Reference Allele Activity')

data_frame(x = seq(-2, 2, by = .01), 
           lik = map_dbl(x, ~likFun(c(.x, 1, 1.5))), 
           prior = map_dbl(x, ~priorFun(c(.x, 1, 1.5))), 
           post = map_dbl(x, ~postFun(c(.x, 1, 1.5)))) %>% 
  gather(dist, dens, -x) %>% 
  mutate(dens = exp(dens + 1000)) %>% #have to add 1000 because exp(-very large number) always evaluates to 0. It gets normalized out later.
  group_by(dist) %>%
  mutate(dens = dens / (sum(dens))) %>% 
  ungroup %>% 
  ggplot(aes(x, dens, color = dist)) + 
  geom_line() + 
  ylab('density') + 
  xlab('Reference Allele Activity') + 
  ggtitle('chr:6 pos:109625879 MPRApos:2/3 G -> A, mutMu = 1, sigma = 1.5')

#load('~/bayesianMPRA/outputs/30NNmcmcOutputs/6_109625879_2-3.RData')
load('~/bayesianMPRA/outputs/30NNmcmcOutputs/22_32875190_1-3.RData')
HDIofMCMC = function(sampleVec, credMass=0.95) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
} # From Kruschke

mcmcres = res$batch %>% 
  as.data.frame() %>% 
  as.tbl %>% 
  mutate(TS = V2 - V1)
hdi = HDIofMCMC(mcmcres$V2 - mcmcres$V1, .99)

mcmcres %>% 
  ggplot(aes(TS)) + 
  geom_histogram(bins = 50) +
  geom_segment(aes(x = hdi[1],
                       y = -5,
                       xend = hdi[2],
                       yend = -5), inherit.aes = FALSE,
               color = 'purple',
               lwd = 1.5) +
  ggtitle('MCMC posterior of Transcriptional Shift \nfor chr22:32875190 G → A') +
  labs(subtitle = expression(italic('with 99% HDI'))) +
  xlab('Transcriptional Shift') + 
  theme(text = element_text(size = 20),
        plot.subtitle = element_text(size = 15))

library(plotly)

mat = expand.grid(seq(0,1, by = .05), seq(.8, 1.7, by = .05)) %>% 
  as.tbl %>% 
  mutate(likDens = exp(map2_dbl(Var1, Var2, ~likFun(c(.x, .y, 1.5))) + 480))

x = seq(0,2, length.out = 40)
y = seq(-1, 3, length.out = 40)
z = matrix(1:(40**2), nrow = 40)

for (i in 1:40) {
  for (j in 1:40) {
    z[i,j] = exp(likFun(c(x[i], y[j], 1.6)))
  }
}
f <- list(
  family = "Courier New, monospace",
  size = 18,
  color = "#7f7f7f"
)
p = plot_ly(x = x,
            y = y,
            z = z,
            type = 'surface') %>% 
  layout(scene = list(xaxis = list(title = 'refMu')))
p

