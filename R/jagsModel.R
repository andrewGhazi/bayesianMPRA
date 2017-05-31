# Could we just define a JAGS model that works on the raw counts using a negative binomial model rather than the activity based stuff
# This would also naturally account for variants with barcodes at count 0

library(arm)
library(rjags)
library(coda)

set.seed(1280)
x = runif(50)
y = 3 + 2*x + rnorm(50, 0,2)

dat = MPRA %>% filter(construct == tmp)


modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ 
  }
}


"

initList = c();

jagsModel = jags.model(file = textConnection(modelString),
                       data = dataList,
                       inits = initList,
                       n.chains = 3, 
                       n.adapt = 500)

dir="/mnt/labhome/andrew/MPRA/paper_data/"

exampleData = read_delim(file = paste0(dir, "Raw/", "RBC_MPRA_minP_raw.txt"),
                         delim = "\t",
                         col_names = T,
                         col_types = cols(chr = "c")) %>% 
  select(matches('construct|CTRL|DNA|type')) %>% 
  filter(construct == '1 155271258 1/3') %>% 
  mutate(bcIndex = 1:nrow(.)) %>% 
  gather(block, count, -construct, -bcIndex, -type)

