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
