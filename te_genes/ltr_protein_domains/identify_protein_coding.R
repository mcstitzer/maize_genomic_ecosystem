library(rtracklayer)
library(plyr)
library(stringr)
library(dplyr)


al=import.gff3('B73V4.pseudomolecule.ltrdigest.ALLsubtracted.MCSnames.20May2017.gff3.gz')

ltr=data.frame(ID=names(table(al$mcs)))

## get info about each TE
al$prot=str_split_fixed(al$name, '_', 2)[,1]

ltr$GAG=ltr$ID %in% al$mcs[al$prot=='GAG' | al$prot=='GAGCOAT']
ltr$AP=ltr$ID %in% al$mcs[al$prot=='AP']
ltr$INT=ltr$ID %in% al$mcs[al$prot=='INT']
ltr$RT=ltr$ID %in% al$mcs[al$prot=='RT']
ltr$RNaseH=ltr$ID %in% al$mcs[al$prot=='RNaseH']
ltr$ENV=ltr$ID %in% al$mcs[al$prot=='ENV']
ltr$CHR=ltr$ID %in% al$mcs[al$prot=='CHR']

ltr$pol=ltr$AP & ltr$INT & ltr$RT & ltr$RNaseH
ltr$auton= ltr$GAG & ltr$pol

ltrprots=data.frame(TEID=ltr$ID, GAG=ltr$GAG, AP=ltr$AP, INT=ltr$INT, RT=ltr$RT, RNaseH=ltr$RNaseH, ENV=ltr$ENV, CHR=ltr$CHR, pol=ltr$pol, auton=ltr$auton)

write.table(ltrprots, 'B73v4_ltr_proteins.txt', quote=F, sep='\t', row.names=F, col.names=T)

ltrprots$fam=substr(ltrprots$TEID,1,8)
ltrprots$sup=substr(ltrprots$TEID,1,3)
ltrprots_fam=ltrprots %>% group_by(fam,sup) %>% dplyr::summarize(nGAG=sum(GAG, na.rm=T), nAP=sum(AP, na.rm=T), nINT=sum(INT, na.rm=T), nRT=sum(RT, na.rm=T), nRNaseH=sum(RNaseH, na.rm=T), nENV=sum(ENV,na.rm=T), nCHR=sum(CHR,na.rm=T), npol=sum(pol, na.rm=T), nauton=sum(auton, na.rm=T))
write.table(ltrprots_fam, 'B73v4_ltr_proteins.fam.txt', quote=F, sep='\t', row.names=F, col.names=T)



ltrprots_sup=ltrprots_fam %>% group_by(sup) %>% dplyr::summarize(nGAG=sum(nGAG, na.rm=T), nAP=sum(nAP, na.rm=T), nINT=sum(nINT, na.rm=T), nRT=sum(nRT, na.rm=T), nRNaseH=sum(nRNaseH, na.rm=T), nENV=sum(nENV,na.rm=T), nCHR=sum(nCHR,na.rm=T), npol=sum(npol, na.rm=T), nauton=sum(nauton, na.rm=T))
write.table(ltrprots_sup, 'B73v4_ltr_proteins.sup.txt', quote=F, sep='\t', row.names=F, col.names=T)
