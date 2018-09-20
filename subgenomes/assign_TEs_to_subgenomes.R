library(dplyr)
library(rtracklayer)

## loads TEFILE, TEDISJOINED, GENOMENAME, GENOME
source('../GenomeInfo.R')

## download subgenomes
subgenomegff='B73v4.subgenome_reconstruction.gff3'
if (!file.exists(subgenomegff)) {
#    setInternet2(TRUE)
    download.file('ftp://ftp.gramene.org/pub/gramene/CURRENT_RELEASE/gff3/zea_mays/subGenome_A.gff3' ,'subGenome_A.gff3',method="auto")
    download.file('ftp://ftp.gramene.org/pub/gramene/CURRENT_RELEASE/gff3/zea_mays/subGenome_B.gff3' ,'subGenome_B.gff3',method="auto")
    system('cat subGenome_A.gff3 subGenome_B.gff3 > B73v4.subgenome_reconstruction.gff3')
}
## find TEs that overlap a subgenome using the disjoined TE file
system(paste0('bedtools intersect -names subgenome -wo -a ', TEDISJOINED ,' -b B73v4.subgenome_reconstruction.gff3 > B73v4_TE_subgenome.bedout.txt'))
a=read.table('B73v4_TE_subgenome.bedout.txt', header=F)

a$TEID=substr(a$V9, 4,24)
a$subgenome=sub(".*.Subgenome=", "", a$V18)
a$subgenome=sub(";Ancestral.*","", a$subgenome)

## to confirm we only have the real parts of this terrible string:
## should only have A and B values
table(a$subgenome)


s=data.frame(TEID=a$TEID, sup=substr(a$TEID,1,3), fam=substr(a$TEID,1,8), subgenome=a$subgenome)

table(s$sup, s$subgenome)/rowSums(table(s$sup, s$subgenome))

write.table(s, paste0(GENOME, "_TE_subgenome.", Sys.Date(), ".txt"), quote=F, sep='\t', col.names=T, row.names=F)

s_fam=s %>% group_by(fam, sup) %>% dplyr::summarize(A=sum(subgenome=='A'), B=sum(subgenome=='B'))
write.table(s_fam, paste0(GENOME, "_TE_subgenome.fam.", Sys.Date(), ".txt"), quote=F, sep='\t', col.names=T, row.names=F)

s_sup=s_fam %>% group_by(sup) %>% dplyr::summarize(A=sum(A, na.rm=T), B=sum(B, na.rm=T))
write.table(s_sup, paste0(GENOME, "_TE_subgenome.sup.", Sys.Date(), ".txt"), quote=F, sep='\t', col.names=T, row.names=F)

### report genome-wide proportions of bp in each subgenome
gwg=import.gff3('B73v4.subgenome_reconstruction.gff3')

gw=sum(width(gwg)[gwg$Subgenome=='A'])/(sum(width(gwg)[gwg$Subgenome=='A'])+sum(width(gwg)[gwg$Subgenome=='B']))
write.table(gw, paste0('genome-wide-proportion-subgenomeA.', Sys.Date(), '.txt'), quote=F, sep='\t', col.names=F, row.names=F)

