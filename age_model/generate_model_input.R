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
techar=fread('../te_characteristics/B73_TE_individual_copies.2019-10-18.txt')
################### genes nearby
gene=fread('../genes/B73_closest_gene.2019-10-21.txt')
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
#tegenes=fread('../te_genes/proteins/B73.hmmprotein.txt')
#ltrgenes=fread('../te_genes/ltr_protein_domains/B73v4_ltr_proteins.txt')
#ltrgenes$TEID=gsub('B73v4', 'Zm00001d', ltrgenes$TEID)
#ind=merge(ind, tegenes, all.x=T, by=c('TEID', 'fam', 'sup'))
#ind=merge(ind, ltrgenes, all.x=T, by=c('TEID'))
#ind[,c('helprot', 'rveprot', 'tpaseprot', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV', 'CHR', 'pol', 'auton')][is.na(ind[,c('helprot', 'rveprot', 'tpaseprot', 'GAG', 'AP', 'INT', 'RT', 'RNaseH', 'ENV', 'CHR', 'pol', 'auton')])]=F

tegenes=fread('../te_genes/B73.2019-01-29.te_proteins.txt')
tegenes=merge(ind[,1:3], tegenes, all=T) ## because the only things in tegenes are tes with a protein!
tegenes[is.na(tegenes)]=F
ind=merge(ind, tegenes, all.x=T, by=c('TEID', 'fam', 'sup'))


###################### te orfs
orfs=fread('../te_genes/orfs/B73v4_allTE_orf.txt')
ind=merge(ind, orfs, all.x=T, by=c('TEID'))


##################### base composition
basecomp=fread('../base_composition/B73_TE_methylatable.2019-10-18.txt')
basecomp.flank=fread('../base_composition/B73_TE_methylatable.flank.2019-10-18.txt')
colnames(basecomp.flank)[4:8]=paste0(colnames(basecomp.flank)[4:8], '_1kbflank')
mnase=fread('../mnase/B73_TEandFlank_mnase.2018-10-23.txt')
ind=merge(ind, basecomp, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, basecomp.flank, all.x=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, mnase, all.x=T, by=c('TEID', 'fam', 'sup'))

##################### diversity ********* need to fix some duplicate columns here! (as in both number and bp?)
diversity=fread('../diversity/B73.segregatingsites.TEandFlank.2018-10-23.txt')
ind=merge(ind, diversity, all.x=T, by=c('TEID', 'fam', 'sup'))
ind$segsites=NULL ## this has a per bp column, don't need both
ind$flank_segsites=NULL


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


################## gene expression

tgem=fread('../genes/B73_closest_gene_expression.maizegdbWalley.2019-10-21.txt')
colnames(tgem)[3]='TEID'
ind=merge(ind, tgem, all.x=T, by=c('TEID', 'fam', 'sup', 'closest', 'closestgene', 'closestgenetype', 'closest.syntenic', 'closestgene.syntenic', 'closestgenetype.syntenic'))



################## syntenic gene expression

tges=fread('../genes/B73_closest_syntenic_gene_expression.maizegdbWalley.2019-10-21.txt')
colnames(tges)[3]='TEID'
ind=merge(ind, tges, all.x=T, by=c('TEID', 'fam', 'sup', 'closest', 'closestgene', 'closestgenetype', 'closest.syntenic', 'closestgene.syntenic', 'closestgenetype.syntenic'))


################# te expression
fe=read.table('../te_expression/walley_mean_expr.2019-01-28.txt', header=T)

ind=merge(ind, fe, all.x=T, by=c('fam'))


################# recombination
rec=fread('../recombination/B73_recombination.2018-10-28.txt')
colnames(rec)[1]='TEID'
ind=merge(ind, rec, all.x=T, by=c('TEID', 'fam', 'sup'))
ind$cm=NULL
ind$bp=NULL


################## subgenome
sg=fread('../subgenomes/B73_TE_subgenome.2019-01-30.txt') ## get assigned to A, B, or NA

ind=merge(ind, sg, all.x=T, by=c('TEID', 'fam', 'sup'))




##remove extra te_bp that comes in, and flank_bp
ind$te_bp=NULL
ind$flank_bp=NULL
ind$geneID=NULL


## add famsize 
ind$famsize=as.numeric(table(ind$fam)[ind$fam])

dim(ind)

## remove wrong ages, of course these will be correlated!!
ind.age=ind[,-c('chr', 'start', 'end', 'strand', 'source', 'type', 'closestgene', 'closestgenetype', 'k2p', 'raw', 'tn93', 'lk2p', 'lraw', 'ltn93', 'gtsim', 'tbl', 'age', 'mya', 'ya', 'tblmya')]
ind.age=cbind(data.frame(age=ind$age), ind.age)
ind.tbl=ind[,-c('chr', 'start', 'end', 'strand', 'source', 'type', 'closestgene', 'closestgenetype', 'k2p', 'raw', 'tn93', 'lk2p', 'lraw', 'ltn93', 'gtsim', 'tbl', 'age', 'mya', 'ya', 'tblmya')]
ind.tbl=cbind(data.frame(tbl=ind$tbl), ind.tbl)

write.table(ind.age, paste0(GENOME, '.LTRAGE.allDescriptors.', Sys.Date(), '.txt'), quote=F, col.names=T, row.names=F, sep='\t')
write.table(ind.tbl, paste0(GENOME, '.TBLAGE.allDescriptors.', Sys.Date(), '.txt'), quote=F, col.names=T, row.names=F, sep='\t')








