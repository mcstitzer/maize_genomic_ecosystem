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

tgem=fread('../genes/B73_closest_gene_expression.maizegdbWalley.txt')
colnames(tgem)[3]='TEID'
ind=merge(ind, tgem, all.x=T, by=c('TEID', 'fam', 'sup', 'closestgene', 'closestgenetype'))


################# te expression
#************************************* THIS WILL BE UPDATED!!!!!! *********************************
fe=read.table('~/projects/b73_ecology/RawData/expression/family_rpm_devatlas_27Feb.txt', header=T)
s=fe
## combine expression values
fe[,'endo_12d_combo']=rowSums(s[,1:12], na.rm=T)
fe[,'peri_aleurone_combo']=rowSums(s[,13:24], na.rm=T)
fe[,'endocrown_combo']=rowSums(s[,25:36], na.rm=T)
fe[,'symdivz_combo']=rowSums(s[,37:48], na.rm=T)
fe[,'stodivz_combo']=rowSums(s[,49:60], na.rm=T)
fe[,'groz_combo']=rowSums(s[,61:72], na.rm=T)
fe[,'emb_20d_combo']=rowSums(s[,73:84], na.rm=T)
fe[,'embryo_combo']=rowSums(s[,85:96], na.rm=T)
fe[,'germk_combo']=rowSums(s[,97:108], na.rm=T) ## okay i am already tired of typing this - bad idea!
fe[,'leaf_8_combo']=rowSums(s[,109:120], na.rm=T)
fe[,'earp_6m_combo']=rowSums(s[,121:132], na.rm=T)
fe[,'earp_2m_combo']=rowSums(s[,133:144], na.rm=T)
fe[,'root_m_combo']=rowSums(s[,145:156], na.rm=T)
fe[,'root_e_combo']=rowSums(s[,157:168], na.rm=T)
fe[,'root_c_combo']=rowSums(s[,169:180], na.rm=T)
fe[,'p_root_combo']=rowSums(s[,181:192], na.rm=T)
fe[,'s_root_combo']=rowSums(s[,193:204], na.rm=T)
fe[,'pollen_combo']=rowSums(s[,205:216], na.rm=T)
fe[,'silk_combo']=rowSums(s[,217:228], na.rm=T)
fe[,'f_spike_combo']=rowSums(s[,229:240], na.rm=T)
fe[,'inode_6_combo']=rowSums(s[,241:252], na.rm=T)
fe[,'inode_7_combo']=rowSums(s[,253:264], na.rm=T)
fe[,'meristem_combo']=rowSums(s[,265:272], na.rm=T) ## there are only two here - did I do something wrong? - Nope, there are only two in the publication
fe=fe[,-c(1:272)]
fe$sup=substr(rownames(fe), 1,3)
fe$fam=rownames(fe)

colnames(fe)[1:23]=paste0('te_', colnames(fe)[1:23])

ind=merge(ind, fe, all.x=T, by=c('fam', 'sup'))


################# recombination
rec=fread('../recombination/B73_recombination.2018-10-28.txt')
colnames(rec)[1]='TEID'
ind=merge(ind, rec, all.x=T, by=c('TEID', 'fam', 'sup'))
ind$cm=NULL
ind$bp=NULL


##remove extra te_bp that comes in, and flank_bp
ind$te_bp=NULL
ind$flank_bp=NULL

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








