library(rtracklayer)
library(data.table)
library(dplyr)

source('../GenomeInfo.R')

te=fread(paste0(GENOME, '.segsites.txt'), sep='\t')
te$TEID=substr(te$V9,4,24)
flank=fread(paste0(GENOME, '.segsites.flank.txt'), sep='\t')
flank$TEID=substr(flank$V9,4,24)

## get TE lengths to divide by
tl=as.data.frame(te) %>% group_by(TEID) %>% summarize(te_bp=sum(V5-V4+1, na.rm=T), segsites=sum(V10, na.rm=T))

tl$segsites.bp=tl$segsites/tl$te_bp

tl$fam=substr(tl$TEID,1,8)
tl$sup=substr(tl$TEID,1,3)

## although flank should be 1kb, at ends of chromosomes and contigs it can be shorter!
tf=as.data.frame(flank) %>% group_by(TEID) %>% summarize(flank_bp=sum(V5-V4+1, na.rm=T), flank_segsites=sum(V10, na.rm=T))
tf$flank_segsites.bp=tf$flank_segsites/tf$flank_bp

segsites=merge(tl, tf, by='TEID', all=T)

## reorder cols
segsites=segsites[,c('TEID', 'fam', 'sup', 'te_bp', 'segsites', 'segsites.bp', 'flank_bp', 'flank_segsites', 'flank_segsites.bp')]


write.table(segsites, paste0(GENOME, '.segregatingsites.TEandFlank.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=F)
