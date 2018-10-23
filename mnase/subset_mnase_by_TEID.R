library(rtracklayer)
library(data.table)
library(dplyr)

source('../GenomeInfo.R')

m=fread('B73_MNaseHS.bedout.txt')
m$TEID=substr(m$V9, 4,24)
flank=fread('B73_MNaseHS.flank.bedout.txt')
flank$TEID=substr(flank$V9,4,24)

mr=filter(m, V10=='RP') %>% group_by(TEID) %>% summarize(te_bp=mean(V5-V4+1, na.rm=T), n_root_hs=n(), root_bp=sum(V14))
ms=filter(m, V10=='AP') %>% group_by(TEID) %>% summarize(te_bp=mean(V5-V4+1, na.rm=T), n_shoot_hs=n(), shoot_bp=sum(V14))
mrf=filter(flank, V10=='RP') %>% group_by(TEID) %>% summarize(flank_bp=mean(V5-V4+1, na.rm=T), flank_n_root_hs=n(), flank_root_bp=sum(V14))
msf=filter(flank, V10=='AP') %>% group_by(TEID) %>% summarize(flank_bp=mean(V5-V4+1, na.rm=T), flank_n_shoot_hs=n(), flank_shoot_bp=sum(V14))

mnase=merge(mr, ms, all=T)
mnase=merge(mnase, mrf, all=T)
mnase=merge(mnase, msf, all=T)


mnase$root_prop=mnase$root_bp/mnase$te_bp
mnase$shoot_prop=mnase$shoot_bp/mnase$te_bp
mnase$flank_root_prop=mnase$flank_root_bp/mnase$flank_bp
mnase$flank_shoot_prop=mnase$flank_shoot_bp/mnase$flank_bp

## for the TEs with no overlaps, want to include them in this file
allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')


mnase=merge(mnase, allte[,'TEID'], all=T)
mnase[is.na(mnase)]=0


mnase$sup=substr(mnase$TEID,1,3)
mnase$fam=substr(mnase$TEID,1,8)

## reorder cols
mnase=mnase[,c('TEID', 'sup', 'fam', 'n_root_hs', 'root_bp', 'root_prop', 'n_shoot_hs', 'shoot_bp', 'shoot_prop', 'flank_n_root_hs', 'flank_root_bp', 'flank_root_prop', 'flank_n_shoot_hs', 'flank_shoot_bp', 'flank_shoot_prop')]
write.table(mnase, paste0(GENOME, '_TEandFlank_mnase.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

