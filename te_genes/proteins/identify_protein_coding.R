library(rtracklayer)
library(plyr)
library(stringr)
library(dplyr)
library(rhmmer)

source('../../GenomeInfo.R')

## read these in with rhmmer parser
a=read_tblout('B73.allteorf.out')
h=read_tblout('B73.helteorf.out')

a$TEID=str_split_fixed(a$query_name, '::', 4)[,2]
a$fam=substr(a$TEID, 1,8)
a$sup=substr(a$TEID, 1,3)
a$pfam=str_split_fixed(a$domain_accession, '\\.', 3)[,1]
#head(a[a$sup %in% c('DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX'),])

h$TEID=str_split_fixed(h$query_name, '::', 4)[,2]
h$fam=substr(h$TEID, 1,8)
h$sup=substr(h$TEID, 1,3)
h$pfam=str_split_fixed(h$domain_accession, '\\.', 3)[,1]
#head(a[a$sup %in% c('DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX'),])


#write.table(tir, 'tir_prot_domain.txt', sep='\t', quote=F, row.names=F, col.names=T)

helprot=h %>% group_by(TEID, sup, fam) %>% summarize(helprot=n()>0)

rveprot=a %>% group_by(TEID, sup, fam) %>% summarize(rveprot=any(pfam=='PF00665')) ## rve is the integrase domain, found in retrotransposons and retroviruses https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5751038

tpaseprot=a %>% group_by(TEID, sup, fam) %>% summarize(tpaseprot=any(pfam!='PF00665'))

prot=merge(helprot, rveprot, all=T)
prot=merge(prot, tpaseprot, all=T)

prot$helprot[is.na(prot$helprot)]=F
prot$rveprot[is.na(prot$rveprot)]=F
prot$tpaseprot[is.na(prot$tpaseprot)]=F

write.table(prot, paste0(GENOME, '.hmmprotein.txt'), col.names=T, row.names=F, sep='\t', quote=F)
