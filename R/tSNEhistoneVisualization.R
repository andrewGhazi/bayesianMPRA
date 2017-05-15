#Now that we've got the histone mark data, let's try to produce a tSNE visualization of the histone marking on the Ulirsch variants

library(tsne)
library(Rtsne)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(data.table)
library(tidyr)
library(dtplyr)
library(pbapply)
library(GenomicRanges)

load('~/MPRA/data/varInfo.RData') #loads varInfo
load('~/Qual/data/histoneMarkFileNames.RData') #loads files, a variable giving the file names of the different histone marks


#### Let's get the scores for one histone mark type at first
H3k4me3 = fread('/mnt/bigData2/andrew/histoneENCODE/wiggles/wgEncodeSydhHistoneK562H3k4me3bUcdSig.wiggle', 
            sep = '\t',
            header = FALSE,
            col.names = c('chr', 'start', 'stop', 'value')) %>% tbl_dt

markRanges = GRanges(seqnames = Rle(H3k4me3$chr), 
              ranges = IRanges(start = H3k4me3$start,
                               end = H3k4me3$stop),
              markValue = H3k4me3$value) %>% 
  split(., seqnames(.))

histoneMarks = files %>% gsub('wgEncode|Histone|StdSig|UcdSig', '', .)

getVariantHistoneValue = function(chrStr, pos){
  chrStr %<>% as.character %>% paste0('chr', .)
  findValue = markRanges[which(start(tmp[chrStr]) < pos & end(tmp[chrStr]) > pos)] %>% .[[1]] %>% mcols %>% .$markValue
  if(length(findValue) == 0){
    return(0)
  } else if (length(findValue) > 1){
    stop(paste0('Multiple values found for ', chrStr, ' at position ', pos))
  } else{
    return(findValue)
  }
}

#This finds the appropriate histone mark values for each allele
for(i in 1:length(histoneMarks)){
  histoneDatI = fread(paste0('/mnt/bigData2/andrew/histoneENCODE/wiggles/', files[i], '.wiggle'), 
                  sep = '\t',
                  header = FALSE,
                  col.names = c('chr', 'start', 'stop', 'value')) %>% tbl_dt
  
  markRanges = GRanges(seqnames = Rle(histoneDatI$chr), 
                       ranges = IRanges(start = histoneDatI$start,
                                        end = histoneDatI$stop),
                       markValue = histoneDatI$value) %>% 
    split(., seqnames(.))
  
  varInfo[[histoneMarks[i]]] = sapply(1:nrow(varInfo), function(x){
    return(getVariantHistoneValue(varInfo$chr[x], varInfo$pos[x]))
  })
}

save(varInfo, file = 'data/varInfoWithHistoneMarkAnnotations.RData')

varInfo %>% 
  select(contains('K562')) %>% 
  gather(mark, value) %>% 
  ggplot(., aes(value)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(~mark)

markMat = varInfo %>% 
  select(contains('K562')) %>% 
  distinct() %>% 
  as.data.frame %>% 
  as.matrix

markMat = varInfo %>% 
  select(contains('K562')) %>% 
  select(contains('Broad')) %>% 
  as.data.frame %>% 
  as.matrix


kVals = seq(1,31, by = 2)
kMeansError = vector('numeric', length(kVals))
for(i in 1:length(kVals)){
  kMeansError[i] = kmeans(markMat, centers = kVals[i], iter.max = 20)$tot.withinss
}
plot(kMeansError) 

# This suggests there are roughly 5 clusters. Given there are 1383 distinct
# rows, this means a perplexity of ~250 I guess?


library(ggfortify)
histoneMark_PCA = prcomp(markMat,
                         center = TRUE,
                         scale. = TRUE)
autoplot(histoneMark_PCA)
# PCA doesn't work too well

par(mfrow = c(3,3))
for(p in seq(10, 330, by = 40)){
  histoneMark_tsne = Rtsne(markMat,
                           perplexity = p,
                           max_iter = 800,
                           check_duplicates = FALSE, verbose = TRUE)
  histoneMark_tsne$Y %>% plot(main = paste0('perplexity = ', p))
}

#Doesn't look amazing but you do pretty consistently get some clusters and a
#central cloud. Not sure what it means but it's something.

histoneDist = dist(markMat)

histoneMark_tsne = Rtsne(markMat,
                         perplexity = p,
                         max_iter = 800,
                         check_duplicates = FALSE, verbose = TRUE)
tsneTbl = histoneMark_tsne$Y %>% 
  as.data.frame %>% 
  as.tbl()

ggplot(tsneTbl, aes(V1, V2)) + 
  geom_point() + 
  ylab('') +
  xlab('') + 
  ggtitle('t-SNE visualization of 11 histone mark annotations for >7000 variants')+
  theme(plot.title = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())
ggsave('~/Qual/tSNE/11histoneMarks.png',
       height = 5.32,
       width = 12)
