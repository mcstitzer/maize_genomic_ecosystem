library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)
source('../GenomeInfo.R')
source('color_palette.R')

##### RESULTS
## Features of the TE section

## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
colnames(gene)[3]='TEID'
ind=merge(techar, gene, all=T)
ind$ingene=ind$closest==0
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!


data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(bp=sum(tebp), meanbp=mean(tebp), medianbp=median(tebp), disrupted=sum(pieces!=1)/n(), disruptor=sum(disruptor>1)/n()))

## nonuniform chr, results paragraph 2
data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(meangenedist=mean(closest, na.rm=T), mediangenedist=median(closest, na.rm=T)))
# for nonLTR as a whole!
data.frame(ind %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(meangenedist=mean(closest, na.rm=T), mediangenedist=median(closest, na.rm=T)))



## gene expression - needs reading in a ton more data!!!

### walley expr atlas from maizegdb
tgem=fread('../genes/B73_closest_gene_expression.maizegdbWalley.txt')
## gene rpkm
data.frame(tgem %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(mean=mean(gene_median, na.rm=T), median=median(gene_median, na.rm=T)))
# for nonLTR as a whole!
data.frame(tgem %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(mean=mean(gene_median, na.rm=T), median=median(gene_median, na.rm=T)))
## gene rpkm, genes within 2kb of TE
data.frame(tgem[tgem$closest<2000,] %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(mean=mean(gene_median, na.rm=T), median=median(gene_median, na.rm=T)))
# for nonLTR as a whole!
data.frame(tgem[tgem$closest<2000,] %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(mean=mean(gene_median, na.rm=T), median=median(gene_median, na.rm=T)))



### paragraph 3 ages

tbl=fread('../te_age/tbl_age/B73_terminalbranchlength2018-10-27.txt')
ltr=fread('../te_age/B73v4_recovered_ages.txt')
ltr$TEID=gsub('B73v4', 'Zm00001d', ltr$tename)
ltr$tename=NULL

ind=merge(ind, tbl, all=T, by=c('TEID', 'fam', 'sup'))
ind=merge(ind, ltr, all=T, by='TEID')
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

## convert to real numbers
ind$age=ind$tbl
ind$age[!is.na(ind$k2p)]=ind$k2p[!is.na(ind$k2p)]
ind$mya=ind$age/3.3e-8/2/1e6
ind$ya=ind$age/3.3e-8/2
ind$tblmya=ind$tbl/3.3e-8/2/1e6

## for ltr retros as a whole
data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(meanmya=mean(mya, na.rm=T), medianmya=median(mya, na.rm=T)))
# for nonLTR as a whole!
data.frame(ind %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(meanmya=mean(mya, na.rm=T), medianmya=median(mya, na.rm=T)))

######## 
## how many families of LTR retro are single copy?
sum(table(ind$fam[substr(ind$sup,1,2)=='RL'])==1)/length(table(ind$fam[substr(ind$sup,1,2)=='RL']))
## what is median LTR retro family size
median(table(ind$fam[substr(ind$sup,1,2)=='RL']))
## what is mean LTR retro family size?
mean(table(ind$fam[substr(ind$sup,1,2)=='RL']))
## what is max LTR retrofamily size?
max(table(ind$fam[substr(ind$sup,1,2)=='RL']))

## how many families of TIR are single copy?
sum(table(ind$fam[substr(ind$sup,1,2)=='DT'])==1)/length(table(ind$fam[substr(ind$sup,1,2)=='DT']))
## what is median TIR family size
median(table(ind$fam[substr(ind$sup,1,2)=='DT']))
## what is mean TIR family size?
mean(table(ind$fam[substr(ind$sup,1,2)=='DT']))
## what is max TIR family size?
max(table(ind$fam[substr(ind$sup,1,2)=='DT']))

## how many families of helitron are single copy?
sum(table(ind$fam[substr(ind$sup,1,2)=='DH'])==1)/length(table(ind$fam[substr(ind$sup,1,2)=='DH']))
## what is median helitron family size
median(table(ind$fam[substr(ind$sup,1,2)=='DH']))
## what is mean helitron family size?
mean(table(ind$fam[substr(ind$sup,1,2)=='DH']))
## what is max helitron family size?
max(table(ind$fam[substr(ind$sup,1,2)=='DH']))

## how many families of nonLTR are single copy?
sum(table(ind$fam[substr(ind$sup,1,2)%in% c('RS', 'RI')])==1)/length(table(ind$fam[substr(ind$sup,1,2)%in% c('RS', 'RI')]))
## what is median nonLTR family size
median(table(ind$fam[substr(ind$sup,1,2)%in% c('RS', 'RI')]))
## what is mean nonLTR family size?
mean(table(ind$fam[substr(ind$sup,1,2)%in% c('RS', 'RI')]))
## what is max nonLTR family size?
max(table(ind$fam[substr(ind$sup,1,2)%in% c('RS', 'RI')]))


###########
## family ages - per copy!
data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(underhalfmil=sum(mya<=0.5, na.rm=T)/n(), underhalfmiltbl=sum(tblmya<=0.5, na.rm=T)/n()))
data.frame(ind %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(underhalfmil=sum(mya<=0.5, na.rm=T)/n(), underhalfmiltbl=sum(tblmya<=0.5, na.rm=T)/n()))

## family ages - per family!
data.frame(ind) %>% group_by(sup, fam) %>%dplyr::summarize(anyhalfmil=any(mya<=0.5, na.rm=T), anyhalfmiltbl=any(tblmya<=0.5, na.rm=T)) %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n(), underhalfmiltbl=sum(anyhalfmiltbl, na.rm=T)/n())
data.frame(ind) %>% group_by(sup, fam)%>%  dplyr::summarize(anyhalfmil=any(mya<=0.5, na.rm=T), anyhalfmiltbl=any(tblmya<=0.5, na.rm=T))  %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n(), underhalfmiltbl=sum(anyhalfmiltbl, na.rm=T)/n())



### coding - make sure you've fixed the protein file so that the autonomous ones are in the age_model produced ind!!!

fp=ind %>% group_by(sup, fam) %>% dplyr::summarize(sum(auton)/n())
ind %>% group_by(sup, fam) %>% group_by(substr(sup,1,2))%>% dplyr::summarize(sum(auton)/n())
ind %>% group_by(sup, fam) %>% group_by(substr(sup,1,2)%in% c('RI', 'RS'))%>% dplyr::summarize(sum(auton)/n())


