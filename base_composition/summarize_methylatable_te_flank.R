library(rtracklayer)
library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)

source('../GenomeInfo.R')
te=TEDISJOINED
me=fread('B73v4.filtered_disjoined_TEs.ecology.cytpatterns.out')
me$TEID=substr(me$"9_usercol", 4,24)

colnames(me)[1:18]=c('chr', 'a', 'b', 'start', 'end', 'c', 'd', 'e', 'id', 'p_at', 'p_gc', 'A', 'C', 'G', 'T', 'N', 'other', 'seqlen')
colnames(me)[19:40]=c('CG', 'CAG', 'CTG', 'CCG', 'CAA', 'CAT', 'CAC', 'CTA', 'CTT', 'CTC', 'CCC', 'CCA', 'CCT', 'TTG', 'ATG', 'GTG', 'TAG', 'AAG', 'GAG', 'GGG', 'TGG', 'AGG')

mm=me %>% group_by(TEID) %>% summarize_if(.predicate=function(x) is.integer(x), .funs=funs('sum'))
mm$start=NULL
mm$end=NULL
mm$sup=substr(mm$TEID, 1,3)
mm$fam=substr(mm$TEID, 1,8)


mm_sup=mm %>% group_by(sup) %>% summarize_if(.predicate=function(x) is.integer(x), .funs=funs('sum'))
mm_fam=mm %>% group_by(fam) %>% summarize_if(.predicate=function(x) is.integer(x), .funs=funs('sum'))
  
### summarize in proportions - GC, CG, CHG, CHH
mm$percGC=(mm$C+mm$G)/(mm$seqlen-mm$N)
mm$nCG=mm$CG/(mm$seqlen-mm$N)
mm$nCHG=(mm$CAG + mm$CTG + mm$CCG)/(mm$seqlen-mm$N)  ## these are normalized, because, for example, a CHG could be double counted in CHG and CG if it were CCG
mm$nCHH=(mm$CAA + mm$CAT + mm$CAC + mm$CTA + mm$CTT + mm$CTC + mm$CCC + mm$CCA + mm$CCT + mm$TTG + mm$ATG + mm$GTG + mm$TAG + mm$AAG + mm$GAG + mm$GGG + mm$TGG + mm$AGG)/(mm$seqlen-mm$N)/2

mm_fam$percGC=(mm_fam$C+mm_fam$G)/(mm_fam$seqlen-mm_fam$N)
mm_fam$nCG=mm_fam$CG/(mm_fam$seqlen-mm_fam$N)
mm_fam$nCHG=(mm_fam$CAG + mm_fam$CTG + mm_fam$CCG)/(mm_fam$seqlen-mm_fam$N)  ## these are normalized, because, for example, a CHG could be double counted in CHG and CG if it were CCG
mm_fam$nCHH=(mm_fam$CAA + mm_fam$CAT + mm_fam$CAC + mm_fam$CTA + mm_fam$CTT + mm_fam$CTC + mm_fam$CCC + mm_fam$CCA + mm_fam$CCT + mm_fam$TTG + mm_fam$ATG + mm_fam$GTG + mm_fam$TAG + mm_fam$AAG + mm_fam$GAG + mm_fam$GGG + mm_fam$TGG + mm_fam$AGG)/(mm_fam$seqlen-mm_fam$N)/2

mm_sup$percGC=(mm_sup$C+mm_sup$G)/(mm_sup$seqlen-mm_sup$N)
mm_sup$nCG=mm_sup$CG/(mm_sup$seqlen-mm_sup$N)
mm_sup$nCHG=(mm_sup$CAG + mm_sup$CTG + mm_sup$CCG)/(mm_sup$seqlen-mm_sup$N)  ## these are normalized, because, for example, a CHG could be double counted in CHG and CG if it were CCG
mm_sup$nCHH=(mm_sup$CAA + mm_sup$CAT + mm_sup$CAC + mm_sup$CTA + mm_sup$CTT + mm_sup$CTC + mm_sup$CCC + mm_sup$CCA + mm_sup$CCT + mm_sup$TTG + mm_sup$ATG + mm_sup$GTG + mm_sup$TAG + mm_sup$AAG + mm_sup$GAG + mm_sup$GGG + mm_sup$TGG + mm_sup$AGG)/(mm_sup$seqlen-mm_sup$N)/2

write.table(mm, 'B73v4_TE_methylatable.txt', quote=F, sep='\t', row.names=F, col.names=T)
write.table(mm_sup, 'B73v4_TE_methylatable.superfam.txt', quote=F, sep='\t', row.names=F, col.names=T)
write.table(mm_fam, 'B73v4_TE_methylatable.fam.txt', quote=F, sep='\t', row.names=F, col.names=T)


######################
## FLANKING ##########
######################

me=fread('B73v4.filtered_disjoined_TEs.flank.cytpatterns.out')
me$TEID=substr(me$"9_usercol", 4,24)

colnames(me)[1:18]=c('chr', 'a', 'b', 'start', 'end', 'c', 'd', 'e', 'id', 'p_at', 'p_gc', 'A', 'C', 'G', 'T', 'N', 'other', 'seqlen')
colnames(me)[19:40]=c('CG', 'CAG', 'CTG', 'CCG', 'CAA', 'CAT', 'CAC', 'CTA', 'CTT', 'CTC', 'CCC', 'CCA', 'CCT', 'TTG', 'ATG', 'GTG', 'TAG', 'AAG', 'GAG', 'GGG', 'TGG', 'AGG')

mm=me %>% group_by(TEID) %>% summarize_if(.predicate=function(x) is.integer(x), .funs=funs('sum'))
mm$start=NULL
mm$end=NULL
mm$sup=substr(mm$TEID, 1,3)
mm$fam=substr(mm$TEID, 1,8)


mm_sup=mm %>% group_by(sup) %>% summarize_if(.predicate=function(x) is.integer(x), .funs=funs('sum'))
mm_fam=mm %>% group_by(fam) %>% summarize_if(.predicate=function(x) is.integer(x), .funs=funs('sum'))
  
### summarize in proportions - GC, CG, CHG, CHH
mm$percGC=(mm$C+mm$G)/(mm$seqlen-mm$N)
mm$nCG=mm$CG/(mm$seqlen-mm$N)
mm$nCHG=(mm$CAG + mm$CTG + mm$CCG)/(mm$seqlen-mm$N)  ## these are normalized, because, for example, a CHG could be double counted in CHG and CG if it were CCG
mm$nCHH=(mm$CAA + mm$CAT + mm$CAC + mm$CTA + mm$CTT + mm$CTC + mm$CCC + mm$CCA + mm$CCT + mm$TTG + mm$ATG + mm$GTG + mm$TAG + mm$AAG + mm$GAG + mm$GGG + mm$TGG + mm$AGG)/(mm$seqlen-mm$N)/2

mm_fam$percGC=(mm_fam$C+mm_fam$G)/(mm_fam$seqlen-mm_fam$N)
mm_fam$nCG=mm_fam$CG/(mm_fam$seqlen-mm_fam$N)
mm_fam$nCHG=(mm_fam$CAG + mm_fam$CTG + mm_fam$CCG)/(mm_fam$seqlen-mm_fam$N)  ## these are normalized, because, for example, a CHG could be double counted in CHG and CG if it were CCG
mm_fam$nCHH=(mm_fam$CAA + mm_fam$CAT + mm_fam$CAC + mm_fam$CTA + mm_fam$CTT + mm_fam$CTC + mm_fam$CCC + mm_fam$CCA + mm_fam$CCT + mm_fam$TTG + mm_fam$ATG + mm_fam$GTG + mm_fam$TAG + mm_fam$AAG + mm_fam$GAG + mm_fam$GGG + mm_fam$TGG + mm_fam$AGG)/(mm_fam$seqlen-mm_fam$N)/2

mm_sup$percGC=(mm_sup$C+mm_sup$G)/(mm_sup$seqlen-mm_sup$N)
mm_sup$nCG=mm_sup$CG/(mm_sup$seqlen-mm_sup$N)
mm_sup$nCHG=(mm_sup$CAG + mm_sup$CTG + mm_sup$CCG)/(mm_sup$seqlen-mm_sup$N)  ## these are normalized, because, for example, a CHG could be double counted in CHG and CG if it were CCG
mm_sup$nCHH=(mm_sup$CAA + mm_sup$CAT + mm_sup$CAC + mm_sup$CTA + mm_sup$CTT + mm_sup$CTC + mm_sup$CCC + mm_sup$CCA + mm_sup$CCT + mm_sup$TTG + mm_sup$ATG + mm_sup$GTG + mm_sup$TAG + mm_sup$AAG + mm_sup$GAG + mm_sup$GGG + mm_sup$TGG + mm_sup$AGG)/(mm_sup$seqlen-mm_sup$N)/2

write.table(mm, 'B73v4_TE_methylatable.flank.txt', quote=F, sep='\t', row.names=F, col.names=T)
write.table(mm_sup, 'B73v4_TE_methylatable.flank.superfam.txt', quote=F, sep='\t', row.names=F, col.names=T)
write.table(mm_fam, 'B73v4_TE_methylatable.flank.fam.txt', quote=F, sep='\t', row.names=F, col.names=T)

