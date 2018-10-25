library(data.table)
library(dplyr)
library(fields)

source('../GenomeInfo.R')

## go through each TE and identify methylation across TE and flank, write to file

for ( TISSUE in c('anther', 'earshoot', 'flag_leaf', 'SAM', 'all3')){
  ORIGTISSUE=TISSUE
  if (ORIGTISSUE =='flag_leaf'){TISSUE='flagleaf'}
  if ( !file.exists(paste0('B73v4.Zm00001d.methylation_spread.', TISSUE, '.txt'))){
    print(TISSUE)
  bed=fread(paste0('TE_methylation_overlap.', ORIGTISSUE, '_Zm00001d.bed'))
  bed$fam=substr(bed$V9, 4, 11)
  bed$TEID=substr(bed$V9,4,24)
  bed$sup=substr(bed$V9,4,6)
  bed$cg=as.numeric(bed$V15)
  bed$chg=as.numeric(bed$V21)
  bed$chh=as.numeric(bed$V27)
  mC=bed %>% group_by(TEID) %>% summarize(avg_cg=mean(cg, na.rm=T), avg_chg=mean(chg, na.rm=T), avg_chh=mean(chh, na.rm=T))
  colnames(mC)[2:4]=paste0(TISSUE, '_', colnames(mC)[2:4])
#  write.table(mC, paste0('B73v4.Zm00001d.TE_methylation.', TISSUE, '.txt'), row.names=F, col.names=T, sep='\t', quote=F)
  bed=fread(paste0('methylationbins_closest40.', ORIGTISSUE, '_Zm00001d.bed'))
  bed$fam=substr(bed$V9, 4, 11)
  bed$TEID=substr(bed$V9,4,24)
  bed$sup=substr(bed$V9,4,6)
  bed$abs.dist=abs(bed$V28)
  bed$cg=as.numeric(bed$V15)
  bed$chg=as.numeric(bed$V21)
  bed$chh=as.numeric(bed$V27)
## make data frame to fill in
#all TEs
  TEIDs=unique(bed$te)
  met.data.all <- data.frame(matrix(NA,nrow=length(TEIDs),ncol=60))
  rownames(met.data.all) <- TEIDs
  colnames(met.data.all) <- paste(rep(c("CG","CHG","CHH"),each=20),rep(seq(from=100,to=2000,by=100),3), TISSUE,sep="_")


bed$window=cut(bed$abs.dist, breaks=seq(from=0, to=2000, by=100))   

    
# Use stats.bin() to fill in averages (very slow as a for loop, there is probably a better way that I do not know)

flank.cg=dcast(bed, TEID+fam+sup~window, value.var='cg', fun.aggregate=mean, na.rm=T)    
flank.chg=dcast(bed, TEID+fam+sup~window, value.var='chg', fun.aggregate=mean, na.rm=T)    
flank.chh=dcast(bed, TEID+fam+sup~window, value.var='chh', fun.aggregate=mean, na.rm=T)   
colnames(flank.cg)[4:23]=paste0(TISSUE, '_flank_cg_', 1:20, '00')
flank.cg$"NA"=NULL
colnames(flank.chg)[4:23]=paste0(TISSUE, '_flank_chg_', 1:20, '00')
flank.chg$"NA"=NULL
colnames(flank.chh)[4:23]=paste0(TISSUE, '_flank_chh_', 1:20, '00')
flank.chh$"NA"=NULL
met.data.all=merge(flank.cg, flank.chg, all=T)
met.data.all=merge(met.data.all, flank.chh, all=T)


allmC=merge(mC, met.data.all)
write.table(allmC, paste0(GENOME, '.TEandFlank_methylation.', TISSUE, '.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t', quote=F)

}
}




