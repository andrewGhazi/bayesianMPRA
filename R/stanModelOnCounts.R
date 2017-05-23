model = "
data{
int<lower=0> bcN ;
int<lower=0> DNAcount[bcN] ;
int block[N] ;
int acid[N] ;
}
parameters{

}
model{
  y ~ neg_binomial_2(mu_block, phi_block)
}
generated quantities{
}
"