library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')


## add in te proteins
tegenes=fread('proteins/B73.hmmprotein.txt')
ltrgenes=fread('ltr_protein_domains/B73v4_ltr_proteins.txt')
ltrgenes$TEID=gsub('B73v4', 'Zm00001d', ltrgenes$TEID)

ind=merge(tegenes, ltrgenes, all.x=T, by=c('TEID'))
ind[,c('helprot', 'rveprot', 'tpaseprot', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV', 'CHR', 'pol', 'auton')][is.na(ind[,c('helprot', 'rveprot', 'tpaseprot', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV', 'CHR', 'pol', 'auton')])]=F


ind=merge(ind, ind%>%group_by(fam) %>% summarize(GAGfam=sum(GAG)>0), by='fam')

ind=merge(ind, ind%>%group_by(fam) %>% summarize(APfam=sum(AP)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(INTfam=sum(INT)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(RTfam=sum(RT)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(RNaseHfam=sum(RNaseH)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(ENVfam=sum(ENV)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(CHRfam=sum(CHR)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(polfam=sum(pol)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(autonfam=sum(auton)>0), by='fam')


ind=merge(ind, ind%>%group_by(fam) %>% summarize(helprotfam=sum(helprot)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(rveprotfam=sum(rveprot)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% summarize(tpaseprotfam=sum(tpaseprot)>0), by='fam')

## need to add write to file at the end
write.table(ind, paste0(GENOME, '.te_proteins.txt'), col.names=T, row.names=F, sep='\t', quote=F)


