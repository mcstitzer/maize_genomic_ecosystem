library(rtracklayer)
library(data.table)
library(stringr)
library(dplyr)
library(plyr)


source('../GenomeInfo.R')

te=import.gff3(TEFILE)
te$sup=substr(te$ID,1,3)
te$fam=substr(te$ID, 1,8)

ted=import.gff3(TEDISJOINED)

ted$sup=substr(ted$ID,1,3)
ted$fam=substr(ted$ID, 1,8)
teds=as.data.frame(ted) %>% group_by(ID) %>% dplyr::summarise(bp=sum(width), pieces=n())
te$bp=mapvalues(te$ID, from=teds$ID, to=teds$bp)
te$bp=as.numeric(as.character(te$bp))
te$pieces=mapvalues(te$ID, from=teds$ID, to=teds$pieces)
te$pieces=as.numeric(as.character(te$pieces))
te$disruptor=countOverlaps(te, type='within', ignore.strand=T)

te$disruptor.samefam=sapply(1:length(te), function(x) countOverlaps(te[x,], te[te$fam==te$fam[x],], ignore.strand=T)) ## I think this is counting any piece any depth in the nest - is this what I want?
te$disruptor.samesup=sapply(1:length(te), function(x) countOverlaps(te[x,], te[te$sup==te$sup[x],], ignore.strand=T))

tt=data.frame(TEID=te$ID, chr=seqnames(te), start=start(te), end=end(te), strand=strand(te), source=te$source, type=te$type, sup=te$sup, fam=te$fam, tebp=te$bp, tespan=width(te), pieces=te$pieces, disruptor=te$disruptor, disruptor.samefam=te$disruptor.samefam, disruptor.samesup=te$disruptor.samesup)
write.table(tt, paste0(GENOME, '_TE_individual_copies.', Sys.Date(), '.txt'), quote=F, sep='\t', col.names=T, row.names=F)
