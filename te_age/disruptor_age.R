library(rtracklayer)
library(dplyr)


## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
colnames(gene)[3]='TEID'
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


a=import.gff3('../B73.2018-09-19.allTE.gff3')
sum(ind$pieces==1 & ind$disruptor==2)

a$test=mapvalues(a$ID, from=ind$TEID, to=(ind$pieces==1 & ind$disruptor==2))
a$age=as.numeric(mapvalues(a$ID, from=ind$TEID, to=ind$age))

at=a[a$test & !is.na(a$test),]

ov=findOverlaps(at, a)

at$selfage=as.numeric(mapvalues(at$ID, from=ind$TEID, to=ind$age))

at$externalage=as.numeric(mapvalues(at$ID





