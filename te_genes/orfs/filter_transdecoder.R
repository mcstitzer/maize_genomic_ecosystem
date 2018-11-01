library(rtracklayer)
library(dplyr)

## read in this transdecoder file (which only contains the LONGEST orf for each)
## note that there is a minimum ORF length of 100 aa, so some TEs won't have an ORF called.
a=import.gff3('B73.eachTE.fa.transdecoder.gff3')
b=a[a$type=='CDS',]

## make a new df that only has orf length
orf=data.frame(TEID=as.character(seqnames(b)), orfAA=width(b)/3)

write.table(orf, 'B73v4_allTE_orf.txt', sep='\t', quote=F, row.names=F, col.names=T)

orf$fam=substr(orf$TEID, 1,8)
orf$sup=substr(orf$TEID,1,3)

orf_fam=orf %>% group_by(sup, fam) %>% dplyr::summarize(avg_orfAA=mean(orfAA, na.rm=T))
write.table(orf_fam, 'B73v4_allTE_orf.fam.txt', sep='\t', quote=F, row.names=F, col.names=T)

orf_sup=orf_fam %>% group_by(sup) %>% dplyr::summarize(avg_orfAA=mean(avg_orfAA, na.rm=T))
write.table(orf_sup, 'B73v4_allTE_orf.sup.txt', sep='\t', quote=F, row.names=F, col.names=T)
