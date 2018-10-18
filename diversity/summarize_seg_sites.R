library(rtracklayer)
library(data.table)
library(dplyr)

source('../GenomeInfo.R')

te=import.gff3($TEDISJOINED)

$GENOME.segsites.txt

m=fread(paste0(GENOME, '.segsites.txt'))
m=m[m$V1 %in% 1:10,] ## get rid of contigs
m$TEID=substr(m$V9, 4,24)

m$segsites=m$V10



## get TE lengths to divide by
tl=as.data.frame(te) %>% group_by(ID) %>% summarize(te_bp=sum(width, na.rm=T))
names(tl)[1]='TEID'

ms=full_join(m, tl)
ms$fam=substr(ms$TEID,1,8)
ms$sup=substr(ms$TEID,1,3)

segsites=ms %>% group_by(TEID) %>% summarize(segsites=sum(segsites), te_bp=sum(te_bp), segsites.bp=sum(segsites)/sum(te_bp))
segsites$sup=substr(segsites$TEID,1,3)
segsites$fam=substr(segsites$TEID, 1,8)






m=fread('B73v4_TE_segsites.FLANK1kbEACH.txt')
m=m[m$V1 %in% 1:10,] ## get rid of contigs
m$TEID=substr(m$V9, 4,24)

m$segsites=m$V10

segsites= m %>% group_by(TEID) %>% summarize(segsites_flank=sum(segsites), denominator=sum(V5-V4), segsites.bp_flank=sum(segsites)/sum(V5-V4))

segsites$sup=substr(segsites$TEID,1,3)
segsites$fam=substr(segsites$TEID, 1,8)

