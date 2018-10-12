library(ape)
library(reshape2)
library(stringr)
library(plyr)


filenames=list.files('aligned_ltrs', pattern='.fa')
ages=data.frame(tename=str_split_fixed(filenames, '\\.', 3)[,1])
ages$k2p=NA
ages$raw=NA
ages$tn93=NA

for (i in 1:length(filenames)){
 ltr=read.FASTA(paste0('aligned_ltrs/',filenames[i]))
 if( !is.null(ltr) ){
 d=dist.dna(ltr, model='K80', gamma=F)
 ages$k2p[i]=d
 raw=dist.dna(ltr, model='raw', gamma=F)
 ages$raw[i]=raw
 tn93=dist.dna(ltr, model='tn93', gamma=F)
 ages$tn93[i]=tn93
  }
}


#te=readRDS('~/projects/arabidopsis/nested_ltr/TAIR10_Chr.all.ltrs.mask.GRanges.RDS')
#te=import.gff3('../SOMETHING')
te=read.table('../B73v4_TE_individual_copies.txt', header=T)
ages$lk2p=100-(ages$k2p*100)
ages$lraw=100-(ages$raw*100)
ages$ltn93=100-(ages$tn93*100)
ages$gtsim=as.numeric(as.character(mapvalues(ages$tename, from=te$TEID, to=te$ltr_similarity, warn_missing=F)))



png('ltrharvest_age_difference.%03d.png')
plot(ages$lraw~ages$gtsim, xlab='LTR harvest age', ylab='RAW mafft age', pch=19, type='p')
plot(ages$lk2p~ages$gtsim, xlab='LTR harvest age', ylab='K2P mafft age', pch=19, type='p')
plot(ages$ltn93~ages$gtsim, xlab='LTR harvest age', ylab='TN93 mafft age', pch=19, type='p')

dev.off()

write.table(ages, 'B73v4_recovered_ages.txt', quote=F, col.names=T, row.names=F, sep='\t')
