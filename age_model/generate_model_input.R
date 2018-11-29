library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')

# TE characteristics
# Genes nearby
# ages
# te proteins
# base composition
# diversity
# methylation

################## characteristics
## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
################### genes nearby
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
colnames(gene)[3]='TEID'
################### ages
tbl=fread('../te_age/tbl_age/B73_terminalbranchlength2018-10-27.txt')
ltr=fread('../te_age/B73v4_recovered_ages.txt')
ltr$TEID=gsub('B73v4', 'Zm00001d', ltr$tename)
ltr$tename=NULL

ind=merge(techar, gene, all=T)
ind$ingene=ind$closest==0
ind=merge(ind, ltr, all.x=T, by='TEID')
ind=merge(ind, tbl, all.x=T, by=c('TEID', 'fam', 'sup'))

ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

ind$age=ind$tbl
ind$age[!is.na(ind$k2p)]=ind$k2p[!is.na(ind$k2p)]
ind$mya=ind$age/3.3e-8/2/1e6
ind$ya=ind$age/3.3e-8/2
ind$tblmya=ind$tbl/3.3e-8/2/1e6

##################### te proteins
## add in te proteins
tegenes=fread('../te_genes/proteins/B73.hmmprotein.txt')
ltrgenes=fread('../te_genes/ltr_protein_domains/B73v4_ltr_proteins.txt')
ltrgenes$TEID=gsub('B73v4', 'Zm00001d', ltrgenes$TEID)
ind=merge(ind, tegenes, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, ltrgenes, all.x=T, by=c('TEID'))
ind[,c('helprot', 'rveprot', 'tpaseprot', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV', 'CHR', 'pol', 'auton')][is.na(ind[,c('helprot', 'rveprot', 'tpaseprot', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV', 'CHR', 'pol', 'auton')])]=F

##################### base composition
basecomp=fread('../base_composition/B73_TE_methylatable.txt')
basecomp.flank=fread('../base_composition/B73_TE_methylatable.flank.txt')
colnames(basecomp.flank)[4:7]=paste0(colnames(basecomp.flank)[4:7], '_1kbflank')
mnase=fread('../mnase/B73_TEandFlank_mnase.2018-10-23.txt')
ind=merge(ind, basecomp, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, basecomp.flank, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, mnase, all.x=T, by=c('TEID', 'fam', 'sup'))

##################### diversity ********* need to fix some duplicate columns here!
diversity=fread('../diversity/B73.segregatingsites.TEandFlank.2018-10-23.txt')
ind=merge(ind, diversity, all.x=T, by=c('TEID', 'fam', 'sup'))


##################### methylation
anther=fread('../methylation/B73.TEandFlank_methylation.anther.2018-10-24.txt')
earshoot=fread('../methylation/B73.TEandFlank_methylation.earshoot.2018-10-24.txt')
flagleaf=fread('../methylation/B73.TEandFlank_methylation.flagleaf.2018-10-25.txt')
all3=fread('../methylation/B73.TEandFlank_methylation.all3.2018-10-25.txt')
sam=fread('../methylation/B73.TEandFlank_methylation.SAM.2018-10-25.txt')
#h3k9me2=fread()
ind=merge(ind, anther, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, earshoot, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, flagleaf, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, all3, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, sam, all.x=T, by=c('TEID', 'fam', 'sup'))



dim(ind)

ind.age=ind[,-c('chr', 'start', 'end', 'strand', 'source', 'type', 'closestgene', 'closestgenetype', 'k2p', 'raw', 'tn93', 'lk2p', 'lraw', 'ltn93', 'gtsim', 'tbl', 'age', 'mya', 'ya', 'tblmya')]
ind.age=cbind(data.frame(age=ind$age), ind.age)
ind.tbl=ind[,-c('chr', 'start', 'end', 'strand', 'source', 'type', 'closestgene', 'closestgenetype', 'k2p', 'raw', 'tn93', 'lk2p', 'lraw', 'ltn93', 'gtsim', 'tbl', 'age', 'mya', 'ya', 'tblmya')]
ind.tbl=cbind(data.frame(tbl=ind$tbl), ind.tbl)

write.table(ind.age, paste0(GENOME, '.LTRAGE.allDescriptors.', Sys.Date(), '.txt'), quote=F, col.names=T, row.names=F, sep='\t')
write.table(ind.tbl, paste0(GENOME, '.TBLAGE.allDescriptors.', Sys.Date(), '.txt'), quote=F, col.names=T, row.names=F, sep='\t')








