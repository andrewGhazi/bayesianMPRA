#Pull stuff from the other directories. Qual, plateletMPRA, and MPRA.

library(magrittr)

file.copy('~/Qual/scripts/', '~/bayesianMPRA/', recursive = TRUE)
file.copy('~/Qual/data/varInfoWithHistoneMarkAnnotations.RData', 'data/') #Most recent version of varInfo, describing the 7803 variants tested in Ulirsch 2016

###
QRmd = list.files('~/Qual/') %>% grep('.Rmd',., value = TRUE)
for(i in 1:length(QRmd)){
  file.copy(paste0('~/Qual/', QRmd[i]), 'RMarkdown/')
}

QR = list.files('~/Qual/') %>% grep('.R$',., value = TRUE)
for(i in 1:length(QR)){
  file.copy(paste0('~/Qual/', QR[i]), 'scripts/')
}

file.copy('~/Qual/UCSCscripts/', 'externalSoftware/', recursive = TRUE)

file.copy('~/MPRA/data/gatheredMPRA.RData', 'data/')
