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
ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-01-30.txt')

## could make these use GENOME to get the file name
#techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
#gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
#colnames(gene)[3]='TEID'
#ind=merge(techar, gene, all=T)
#ind$ingene=ind$closest==0
#ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!


data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(bp=sum(tebp), meanbp=mean(tebp), medianbp=median(tebp), medianclosest=median(closest, na.rm=T), disrupted=sum(pieces!=1)/n(), disruptor=sum(disruptor>1)/n()))
data.frame(ind %>% group_by(substr(sup,1,2)%in% c('RI', 'RS')) %>% dplyr::summarize(bp=sum(tebp), meanbp=mean(tebp), medianclosest=median(closest, na.rm=T), medianbp=median(tebp), disrupted=sum(pieces!=1)/n(), disruptor=sum(disruptor>1)/n()))
## na.rm for closest is fine because these TEs can be on contigs, some of which lack genes.

## prop within transcript, age
data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(bp=sum(tebp),ingene=sum(closest==0, na.rm=T)/n(), mya=median(age)/2/3.3e-8/1e6))
data.frame(ind %>% group_by(substr(sup,1,2)%in% c('RI', 'RS')) %>% dplyr::summarize(bp=sum(tebp),ingene=sum(closest==0, na.rm=T)/n(), mya=median(age)/2/3.3e-8/1e6))



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
## gene rpkm, genes within 1kb of TE
data.frame(tgem[tgem$closest<1000,] %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(mean=mean(gene_median, na.rm=T), median=median(gene_median, na.rm=T)))
# for nonLTR as a whole!
data.frame(tgem[tgem$closest<1000,] %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(mean=mean(gene_median, na.rm=T), median=median(gene_median, na.rm=T)))


## and TE expression from sarah
ind$exprMedian=apply(ind[,386:408],1,median, na.rm=T)
ind$exprMedianPerCopy=ind$exprMedian/ind$famsize
ind$exprMedianPerBp=ind$exprMedianPerCopy/ind$te_bp
# and tau
tau=function(x){
    t=sum(1-x/max(x))/(length(x)-1)
  }
ind$te_tau=apply(ind[,386:408], 1, tau)
## te rpkm
data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(mean=mean(exprMedianPerCopy, na.rm=T), median=median(exprMedianPerCopy, na.rm=T)))
# for nonLTR as a whole!
data.frame(ind %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(mean=mean(exprMedianPerCopy, na.rm=T), median=median(exprMedianPerCopy, na.rm=T)))



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
ind$mya=ind$age/3.3e-8/2/1e6
## family ages - per copy!
data.frame(ind %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(underhalfmil=sum(mya<=0.5, na.rm=T)/n(), underhalfmiltbl=sum(tblmya<=0.5, na.rm=T)/n()))
data.frame(ind %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(underhalfmil=sum(mya<=0.5, na.rm=T)/n(), underhalfmiltbl=sum(tblmya<=0.5, na.rm=T)/n()))

## family ages - per family!
data.frame(ind) %>% group_by(sup, fam) %>%dplyr::summarize(anyhalfmil=any(mya<=0.5, na.rm=T), anyhalfmiltbl=any(tblmya<=0.5, na.rm=T)) %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n(), underhalfmiltbl=sum(anyhalfmiltbl, na.rm=T)/n())
data.frame(ind) %>% group_by(sup, fam)%>%  dplyr::summarize(anyhalfmil=any(mya<=0.5, na.rm=T), anyhalfmiltbl=any(tblmya<=0.5, na.rm=T))  %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n(), underhalfmiltbl=sum(anyhalfmiltbl, na.rm=T)/n())

data.frame(ind) %>% group_by(sup, fam) %>% dplyr::summarize(anyhalfmil=any(mya<=0.5, na.rm=T)) %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n())
data.frame(ind) %>% group_by(sup, fam)%>%  dplyr::summarize(anyhalfmil=any(mya<=0.5, na.rm=T)) %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n())

data.frame(ind) %>% group_by(sup, fam) %>% dplyr::summarize(anyhalfmil=any(mya<=0.1, na.rm=T)) %>% group_by(substr(sup,1,2)) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n())
data.frame(ind) %>% group_by(sup, fam)%>%  dplyr::summarize(anyhalfmil=any(mya<=0.1, na.rm=T)) %>% group_by(substr(sup,1,2) %in% c('RI', 'RS')) %>% dplyr::summarize(underhalfmil=sum(anyhalfmil, na.rm=T)/n())



### coding - make sure you've fixed the protein file so that the autonomous ones are in the age_model produced ind!!!

ind$auton[ind$sup=='DHH']=ind$helprot[ind$sup=='DHH'] ## make sure autonomous means what I want it to!
ind$auton[substr(ind$sup,1,2)=='DT']=ind$tpaseprot[substr(ind$sup,1,2)=='DT']
ind$auton[ind$sup %in% c('RIT', 'RIL', 'RST')]=ind$rveprot[ind$sup %in% c('RIT', 'RIL', 'RST')]

ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(GAGfam=sum(GAG)>0), by='fam')

ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(APfam=sum(AP)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(INTfam=sum(INT)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(RTfam=sum(RT)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(RNaseHfam=sum(RNaseH)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(ENVfam=sum(ENV)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(CHRfam=sum(CHR)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(polfam=sum(pol)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(autonfam=sum(auton)>0), by='fam')


ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(helprotfam=sum(helprot)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(rveprotfam=sum(rveprot)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(tpaseprotfam=sum(tpaseprot)>0), by='fam')



fp=ind %>% group_by(sup, fam) %>% dplyr::summarize(propAuton=sum(auton)/n(), famsize=n(), medianORF=median(orfAA, na.rm=T), meanORF=mean(orfAA, na.rm=T))
ind %>% group_by(sup, fam) %>% group_by(substr(sup,1,2))%>% dplyr::summarize(sum(auton)/n())
ind %>% group_by(sup, fam) %>% group_by(substr(sup,1,2)%in% c('RI', 'RS'))%>% dplyr::summarize(sum(auton)/n())

fpltr=ind %>% filter(sup %in% c('RLC', 'RLG', 'RLX')) %>% group_by(sup, fam) %>% dplyr::summarize(propGAG=sum(GAG)/n(), propPol=sum(pol)/n(), propAuton=sum(auton)/n(), famsize=n())

library(stargazer)
stargazer(fp[fp$propAuton==0 & fp$famsize>=10,], summary=F, rownames=F, align=T) ## this makes supp table nonautonomous
stargazer(fp[fp$propAuton>=0.75 & fp$famsize>=10,], summary=F, rownames=F, align=T) ## this makes supp table autonomous


## te expression tables
fe=read.table('../te_expression/walley_mean_expr.2019-01-28.txt', header=T)
fe$famsize=mapvalues(as.character(fe$fam), from=ind$fam, to=ind$famsize)
fe$famsize=as.numeric(fe$famsize)

fe$TEfamMaxPerCopy=apply(fe[,-c(1,25:29)], 1, function(x) max(x, na.rm=T))/fe$famsize
## ask which fams are tissue specific
fe[fe$TEfam_tau>0.99 & fe$famsize>10 & !is.na(fe$TEfam_tau) & fe$TEfamMaxPerCopy>1,]


#### max and min base composition
max(ind.fam[ind.fam$famsize>10,]$percGC)
which.max(ind.fam[ind.fam$famsize>10,]$percGC)
ind.fam[ind.fam$famsize>10,][439,]
median(ind.fam[ind.fam$famsize>10,]$percGC)
min(ind.fam[ind.fam$famsize>10,]$percGC)
which.min(ind.fam[ind.fam$famsize>10,]$percGC)
ind.fam[ind.fam$famsize>10,][639,]
min(ind.fam[ind.fam$famsize>10,]$nCG)
which.min(ind.fam[ind.fam$famsize>10,]$nCG)
which.max(ind.fam[ind.fam$famsize>10,]$nCG)
max(ind.fam[ind.fam$famsize>10,]$nCG)
ind.fam[ind.fam$famsize>10,]$fam[213]
ind.fam[ind.fam$famsize>10,]$fam[642]
min(ind.fam[ind.fam$famsize>10,]$nCHG)
which.min(ind.fam[ind.fam$famsize>10,]$nCHG)
ind.fam[ind.fam$famsize>10,]$fam[215]
which.min(ind.fam[ind.fam$famsize>20,]$nCHG)



############### methylation
                
                         
## just hte numbers across all TEs
mean((ind$SAM_avg_cg + ind$all3_avg_cg + ind$flagleaf_avg_cg + ind$earshoot_avg_cg + ind$anther_avg_cg)/5, na.rm=T)
mean((ind$SAM_avg_chg + ind$all3_avg_chg + ind$flagleaf_avg_chg + ind$earshoot_avg_chg + ind$anther_avg_chg)/5, na.rm=T)
mean((ind$SAM_avg_chh + ind$all3_avg_chh + ind$flagleaf_avg_chh + ind$earshoot_avg_chh + ind$anther_avg_chh)/5, na.rm=T)
## table of tissue and mean methylation
mea=ind %>% group_by(sup) %>% dplyr::summarize(SAM_avg_cg=mean(SAM_avg_cg, na.rm=T), SAM_avg_chg=mean(SAM_avg_chg, na.rm=T), SAM_avg_chh=mean(SAM_avg_chh, na.rm=T), 
                                               all3_avg_cg=mean(all3_avg_cg, na.rm=T), all3_avg_chg=mean(all3_avg_cg, na.rm=T), all3_avg_chh=mean(all3_avg_chh, na.rm=T), 
                                               flagleaf_avg_cg=mean(flagleaf_avg_cg, na.rm=T), flagleaf_avg_chg=mean(flagleaf_avg_chg, na.rm=T), flagleaf_avg_chh=mean(flagleaf_avg_chh, na.rm=T),
                                              earshoot_avg_cg=mean(earshoot_avg_cg, na.rm=T), earshoot_avg_chg=mean(earshoot_avg_chg, na.rm=T), earshoot_avg_chh=mean(earshoot_avg_chh, na.rm=T),
                                              anther_avg_cg=mean(anther_avg_cg, na.rm=T), anther_avg_chg=mean(anther_avg_chg, na.rm=T), anther_avg_chh=mean(anther_avg_chh, na.rm=T))

meas=mea %>% group_by(sup) %>% dplyr::summarize(avg_cg=(SAM_avg_cg + all3_avg_cg + flagleaf_avg_cg + earshoot_avg_cg + anther_avg_cg)/5,
                                                avg_chg=(SAM_avg_chg + all3_avg_chg + flagleaf_avg_chg + earshoot_avg_chg + anther_avg_chg)/5,
                                                avg_chh=(SAM_avg_chh + all3_avg_chh + flagleaf_avg_chh + earshoot_avg_chh + anther_avg_chh)/5)
stargazer(merge(meas, mea), summary=F, rownames=F, align=T) ## this makes supp table methylation

med=ind %>% group_by(sup,fam, famsize) %>% dplyr::summarize(SAM_avg_cg=median(SAM_avg_cg, na.rm=T), SAM_avg_chg=median(SAM_avg_chg, na.rm=T), SAM_avg_chh=median(SAM_avg_chh, na.rm=T), 
                                               all3_avg_cg=median(all3_avg_cg, na.rm=T), all3_avg_chg=median(all3_avg_cg, na.rm=T), all3_avg_chh=median(all3_avg_chh, na.rm=T), 
                                               flagleaf_avg_cg=median(flagleaf_avg_cg, na.rm=T), flagleaf_avg_chg=median(flagleaf_avg_chg, na.rm=T), flagleaf_avg_chh=median(flagleaf_avg_chh, na.rm=T),
                                              earshoot_avg_cg=median(earshoot_avg_cg, na.rm=T), earshoot_avg_chg=median(earshoot_avg_chg, na.rm=T), earshoot_avg_chh=median(earshoot_avg_chh, na.rm=T),
                                              anther_avg_cg=median(anther_avg_cg, na.rm=T), anther_avg_chg=median(anther_avg_chg, na.rm=T), anther_avg_chh=median(anther_avg_chh, na.rm=T))
meds=med %>% group_by(sup,fam, famsize) %>% dplyr::summarize(avg_cg=(SAM_avg_cg + all3_avg_cg + flagleaf_avg_cg + earshoot_avg_cg + anther_avg_cg)/5,
                                                avg_chg=(SAM_avg_chg + all3_avg_chg + flagleaf_avg_chg + earshoot_avg_chg + anther_avg_chg)/5,
                                                avg_chh=(SAM_avg_chh + all3_avg_chh + flagleaf_avg_chh + earshoot_avg_chh + anther_avg_chh)/5)

                         
                         
