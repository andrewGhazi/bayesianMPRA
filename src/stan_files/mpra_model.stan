
data{
  int<lower=0> nRefBarcode ; // number of barcodes in ref allele
  int<lower=0> nMutBarcode ;
  int<lower=0> nDNAblocks ; // number of DNA replicates / blocks / samples / transfections
  int<lower=0> nRNAblocks ;
  int<lower=0> refDNAmat[nRefBarcode, nDNAblocks] ; // MPRA count matrix. Rows = barcodes, columns = samples/blocks
  int<lower=0> refRNAmat[nRefBarcode, nRNAblocks] ;
  int<lower=0> mutDNAmat[nMutBarcode, nDNAblocks] ;
  int<lower=0> mutRNAmat[nMutBarcode, nRNAblocks] ;
  real<lower=0> muRefRNAHyperParams[2, nRNAblocks] ; // gamma hyper-parameters on negative binomial parameters
  real<lower=0> phiRefRNAHyperParams[2, nRNAblocks] ;
  real<lower=0> muMutRNAHyperParams[2, nRNAblocks] ;
  real<lower=0> phiMutRNAHyperParams[2, nRNAblocks] ;
  real<lower=0> muDNAHyperParams[2, nDNAblocks] ;
  real<lower=0> phiDNAHyperParams[2, nDNAblocks] ;
}
parameters{
  real<lower=0> muRefDNA[nDNAblocks] ; //mean parameters for each block in each allele for each nucleic acid
  real<lower=0> muRefRNA[nRNAblocks] ;
  real<lower=0> muMutDNA[nDNAblocks] ;
  real<lower=0> muMutRNA[nRNAblocks] ;
  real<lower=0> phiRefDNA[nDNAblocks] ; // size parameters
  real<lower=0> phiRefRNA[nRNAblocks] ;
  real<lower=0> phiMutDNA[nDNAblocks] ;
  real<lower=0> phiMutRNA[nRNAblocks] ;
}
model{


  for (i in 1:nDNAblocks){

    // negative binomial parameters come from gamma hyper-priors
    muRefDNA[i] ~ gamma(muDNAHyperParams[1, i], muDNAHyperParams[2, i]) ;
    phiRefDNA[i] ~ gamma(phiDNAHyperParams[1, i], phiDNAHyperParams[2, i]) ;
    muMutDNA[i] ~ gamma(muDNAHyperParams[1, i], muDNAHyperParams[2, i]) ; 
    phiMutDNA[i] ~ gamma(phiDNAHyperParams[1, i], phiDNAHyperParams[2, i]) ;

    // count data comes from the specified negative binomial
    refDNAmat[,i] ~ neg_binomial_2(muRefDNA[i], phiRefDNA[i]) ;
    mutDNAmat[,i] ~ neg_binomial_2(muMutDNA[i], phiMutDNA[i]) ;
  }

  for (i in 1:nRNAblocks){
    // negative binomial parameters come from gamma hyper-priors
    muRefRNA[i] ~ gamma(muRefRNAHyperParams[1, i], muRefRNAHyperParams[2, i]) ;
    muMutRNA[i] ~ gamma(muMutRNAHyperParams[1, i], muMutRNAHyperParams[2, i]) ;

    phiRefRNA[i] ~ gamma(phiRefRNAHyperParams[1, i], phiRefRNAHyperParams[2, i]) ;
    phiMutRNA[i] ~ gamma(phiMutRNAHyperParams[1, i], phiMutRNAHyperParams[2, i]) ;

    // count data comes from the specified negative binomial
    refRNAmat[,i] ~ neg_binomial_2(muRefRNA[i], phiRefRNA[i]) ;
    mutRNAmat[,i] ~ neg_binomial_2(muMutRNA[i], phiMutRNA[i]) ;
  }
}

