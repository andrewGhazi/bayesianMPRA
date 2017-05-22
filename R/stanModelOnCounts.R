model = "
data{
int<lower=0> N ;
int<lower=0> count[N] ;
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