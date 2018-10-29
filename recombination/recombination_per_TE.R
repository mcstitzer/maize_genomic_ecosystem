library(MonoPoly)
library(rtracklayer)
library(dtplyr)
library(dplyr)
source('tile-methods.R')

source('../GenomeInfo.R')

te=import.gff3(TEFILE)


map=loadGeneticMap('ogut_fifthcM_map_agpv4.txt')
names(map)=c('marker', 'othername', 'seqnames', 'position', 'cumulative')
## clean up some weirdness
map=map[!is.na(map$position),]
map=map[map$position>0,]
starts=map$position[1:(nrow(map)-1)]
starts.chr=map$seqnames[1:(nrow(map)-1)]
ends=map$position[2:nrow(map)]
ends.chr=map$seqnames[2:nrow(map)]
while(sum(ends<starts & ends.chr==starts.chr)){
map=map[-which(ends<starts & ends.chr==starts.chr)-1,]
starts=map$position[1:(nrow(map)-1)]
starts.chr=map$seqnames[1:(nrow(map)-1)]
ends=map$position[2:nrow(map)]
ends.chr=map$seqnames[2:nrow(map)]
}


while(length(rle(map$seqnames)$lengths)>10){
cs=cumsum(rle(map$seqnames)$lengths)
map=map[-cs[rle(map$seqnames)$lengths<5],]
}



### ready to go
approx=approxGeneticMap(map)
s=predictGeneticMap(approx, te)
ete=te
start(ete)=end(te)
e=predictGeneticMap(approx, ete)

d=e$smoothed_cumulative-s$smoothed_cumulative

## now put this in (genetic maps only have chromosomes, not contigs!)

te$cm=NA
te$cm[as.logical(seqnames(te)%in%1:10)]=d

rec= data.frame(te) %>% group_by(ID) %>% summarize(cm=sum(cm, na.rm=T), bp=sum(width, na.rm=T))
rec$cmmb=rec$cm/rec$bp*1000000


rec$sup=substr(rec$ID, 1,3)
rec$fam=substr(rec$ID,1,8)

write.table(rec, paste0(GENOME, '_recombination.', Sys.Date(), '.txt'), quote=F, row.names=F, col.names=T)
