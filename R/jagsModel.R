library(arm)
library(rjags)
library(coda)

set.seed(1280)
x = runif(50)
y = 3 + 2*x + rnorm(50, 0,2)

dat = data.frame(y,x)

modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ 
  }
}


"

initList = ;

jagsModel = jags.model(file = textConnection(modelString),
                       data = dataList,
                       inits = initList,
                       n.chains = 3, 
                       n.adapt = 500)