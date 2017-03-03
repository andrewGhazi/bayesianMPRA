# Download from UCSC like so:
# wget ftp://hgdownload.cse.ucsc.edu/gbdb/hg38/bbi/wgEncodeReg/wgEncodeRegMarkH3k27ac/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig /mnt/bigData2/andrew/histoneENCODE/

# Convert from bigWig to wiggle like so:
# andrew@smirnov:~/Qual/UCSCscripts$ ~/Qual/UCSCscripts/bigWigToWig /mnt/bigData2/andrew/histonENCODE/bigWigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig ~/Qual/data/wiggles/wgEncodeBroadHistoneK562H3k27acStdSig.wiggle

# Remove the hashmark lines from the wiggle like so
#grep -v '^#' ~/Qual/data/wiggles/wgEncodeBroadHistoneK562H3k27acStdSig.wiggle > ~/Qual/data/wigglesNoHashes/wgEncodeBroadHistoneK562H3k27acStdSigNoHashes.wiggle

library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)

# k562H3k27 = read_tsv('~/Qual/data/wigglesNoHashes/wgEncodeBroadHistoneK562H3k27acStdSigNoHashes.wiggle',
#                      col_names = c('chr', 'start', 'stop', 'value'))
# 
# testRegion = k562H3k27 %>% 
#   filter(chr == 'chr7',
#          start > 80750000,     #A region that's ~200kb downstream of CD36, just picked it out because it had a couple peaks
#          stop < 80900000) %>% 
#   mutate(mids = (start + stop)/2,
#          widths = (stop - start))
# 
# ggplot(testRegion, aes(mids, value)) + 
#   geom_col(width = testRegion$widths)

####Okay this lines up with the data visible on the UCSC browser so let's just download every histone mark dataset
fileNames = system('grep K562H[0-9]k[0-9] /mnt/bigData2/andrew/histoneENCODE/Index\\ of\\ _goldenPath_hg19_database.html', intern = TRUE) %>% 
  grep('.txt.gz', ., value = TRUE) %>% 
  grep('Sig', ., value = TRUE)

fileStarts = fileNames %>% str_locate(., 'wg') %>% .[,1]
fileStops = (fileNames %>% str_locate(., '.txt.gz') %>% .[,1]) - 1
files = fileNames %>% str_sub(., start = fileStarts, end = fileStops)
save(files, file = '~/Qual/data/histoneMarkFileNames.RData')

#download each bigWig
setwd('/mnt/bigData2/andrew/histoneENCODE/bigWigs/')
for(i in 1:length(files)){
  download.file(paste0('ftp://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/', files[i], '.bigWig'), destfile = paste0(files[i], '.bigWig'), quiet = TRUE)
}

# bigWigToWig each file
setwd('/mnt/bigData2/andrew/histoneENCODE/wiggles/')
for(x in 1:length(files)){
  system(paste0('~/Qual/UCSCscripts/bigWigToWig ',
                '/mnt/bigData2/andrew/histoneENCODE/bigWigs/', files[x], '.bigWig ',
                files[x], '.wiggle'))
}

#remove hashmark lines from each wiggle
for(x in 1:length(files)){
  system(paste0('grep -v ^# ',
                files[x], '.wiggle > temp && mv temp ', 
                files[x], '.wiggle'))
}