library(rtracklayer)
library(data.table)
library(dplyr)

source('../GenomeInfo.R')

te=import.gff3(TEDISJOINED)
m=fread('B73v4_MNaseHS.bedout.txt')
m$TEID=substr(m$V9, 4,24)
flank=fread('B73v4_MNaseHS.flank.bedout.txt')
flank$TEID=substr(flank$V9,4,24)

mr=filter(m, V10=='RP') %>% group_by(TEID) %>% summarize(n_root_hs=n(), root_bp=sum(V14))
ms=filter(m, V10=='AP') %>% group_by(TEID) %>% summarize(n_shoot_hs=n(), shoot_bp=sum(V14))

mr$sup=substr(mr$TEID,1,3)
mr$fam=substr(mr$TEID, 1,8)
mr_sup=mr %>% group_by(sup) %>% summarize(avg_root_hs=mean(n_root_hs, na.rm=T), tot_root_hs=sum(n_root_hs, na.rm=T), avg_root_bp=mean(root_bp, na.rm=T), tot_root_bp=sum(root_bp, na.rm=T))
mr_fam=mr %>% group_by(fam) %>% summarize(avg_root_hs=mean(n_root_hs, na.rm=T), tot_root_hs=sum(n_root_hs, na.rm=T), avg_root_bp=mean(root_bp, na.rm=T), tot_root_bp=sum(root_bp, na.rm=T)) %>% arrange(avg_root_hs)

ms$sup=substr(ms$TEID,1,3)
ms$fam=substr(ms$TEID, 1,8)
ms_sup=ms %>% group_by(sup) %>% summarize(avg_shoot_hs=mean(n_shoot_hs, na.rm=T), tot_shoot_hs=sum(n_shoot_hs, na.rm=T), avg_shoot_bp=mean(shoot_bp, na.rm=T), tot_shoot_bp=sum(shoot_bp, na.rm=T))
ms_fam=ms %>% group_by(fam) %>% summarize(avg_shoot_hs=mean(n_shoot_hs, na.rm=T), tot_shoot_hs=sum(n_shoot_hs, na.rm=T), avg_shoot_bp=mean(shoot_bp, na.rm=T), tot_shoot_bp=sum(shoot_bp, na.rm=T)) %>% arrange(avg_shoot_hs)

## need te lengths
tl=as.data.frame(te) %>% group_by(ID) %>% summarize(te_bp=sum(width, na.rm=T))
names(tl)[1]='TEID'

tl_sup=as.data.frame(te) %>% group_by(substr(ID,1,3)) %>% summarize(te_bp=sum(width, na.rm=T))
names(tl_sup)[1]='sup'

tl_fam=as.data.frame(te) %>% group_by(substr(ID,1,8)) %>% summarize(te_bp=sum(width, na.rm=T))
names(tl_fam)[1]='fam'


mnase=full_join(mr, ms)
mnase=full_join(mnase, tl)
mnase$root_prop=mnase$root_bp/mnase$te_bp
mnase$shoot_prop=mnase$shoot_bp/mnase$te_bp
mnase$n_root_hs[is.na(mnase$n_root_hs)]=0
mnase$root_bp[is.na(mnase$root_bp)]=0
mnase$n_shoot_hs[is.na(mnase$n_shoot_hs)]=0
mnase$shoot_bp[is.na(mnase$shoot_bp)]=0
mnase$root_prop[is.na(mnase$root_prop)]=0
mnase$shoot_prop[is.na(mnase$shoot_prop)]=0
mnase$sup=substr(mnase$TEID,1,3)
mnase$fam=substr(mnase$TEID,1,8)

write.table(mnase, 'B73v4_TE_mnase.txt', quote=F, sep='\t', row.names=F, col.names=T)


## these may need fixing to deal with the join and 0 values!!
mnase_sup=full_join(mr_sup, ms_sup)
mnase_sup=full_join(mnase_sup, tl_sup)
mnase_sup$root_prop=mnase_sup$tot_root_bp/mnase_sup$te_bp
mnase_sup$shoot_prop=mnase_sup$tot_shoot_bp/mnase_sup$te_bp
write.table(mnase_sup, 'B73v4_TE_mnase.superfam.txt', quote=F, sep='\t', row.names=F, col.names=T)

mnase_fam=full_join(mr_fam, ms_fam)
mnase_fam=full_join(mnase_fam, tl_fam)
mnase_fam$root_prop=mnase_fam$tot_root_bp/mnase_fam$te_bp
mnase_fam$shoot_prop=mnase_fam$tot_shoot_bp/mnase_fam$te_bp
write.table(mnase_fam, 'B73v4_TE_mnase.fam.txt', quote=F, sep='\t', row.names=F, col.names=T)


##############################
##### flanking stuff
##############################

mfr=filter(flank, V10=='RP') %>% group_by(TEID) %>% summarize(n_root_hs=n(), root_bp=sum(V14))
mfs=filter(flank, V10=='AP') %>% group_by(TEID) %>% summarize(n_shoot_hs=n(), shoot_bp=sum(V14))

mfr$sup=substr(mfr$TEID,1,3)
mfr$fam=substr(mfr$TEID, 1,8)

mnase=full_join(mfr, mfs)

mnase$n_root_hs[is.na(mnase$n_root_hs)]=0
mnase$root_bp[is.na(mnase$root_bp)]=0
mnase$n_shoot_hs[is.na(mnase$n_shoot_hs)]=0
mnase$shoot_bp[is.na(mnase$shoot_bp)]=0

mnase$root_prop=mnase$root_bp/1998 ## each TE has two sides, so this summarizes the 1kb on either side
mnase$shoot_prop=mnase$shoot_bp/1998

mnase$sup=substr(mnase$TEID,1,3)
mnase$fam=substr(mnase$TEID,1,8)
mnase=mnase[,c(1,4,5,2,3,6:9)]
colnames(mnase)[4:ncol(mnase)]=paste(colnames(mnase)[4:ncol(mnase)], 'flank', sep='_')


## add in zeros for all TEs without mnase in flanks
mnaseb=full_join(mnase, data.frame(TEID=tl$TEID))
mnaseb$sup=substr(mnaseb$TEID,1,3)
mnaseb$fam=substr(mnaseb$TEID,1,8)
mnaseb[is.na(mnaseb)]=0


write.table(mnaseb, 'B73v4_TE_mnase.flank.txt', quote=F, sep='\t', row.names=F, col.names=T)
