library(RColorBrewer)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')
source('color_palette.R')

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
basecomp=fread('../base_composition/B73_TE_methylatable.txt')
basecomp.flank=fread('../base_composition/B73_TE_methylatable.flank.txt')
colnames(basecomp.flank)[4:7]=paste0(colnames(basecomp.flank)[4:7], '_1kbflank')
mnase=fread('../mnase/B73_TEandFlank_mnase.2018-10-23.txt')

diversity=fread('../diversity/B73.segregatingsites.TEandFlank.2018-10-23.txt')

## always bring this one in!!!
allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
ind=merge(basecomp, basecomp.flank, all=T)
ind=merge(ind, diversity, all=T)
ind=merge(ind, mnase, all=T)
ind=merge(ind, allte, all=T)
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

nrow(ind)
nrow(basecomp)
nrow(basecomp.flank)

genomewide=fread('../base_composition/B73.basecomp_genomewide.txt')
genomewide.ss=fread('../diversity/genomewide_segsites.txt')
colnames(genomewide.ss)[1]='genomewide_bp_segsites'
genomewide.mnase=read.table(paste0('../mnase/', GENOME, '.genomewide.mnase.txt'))
colnames(genomewide.mnase)=c('n_root_hs_genomewide', 'V2', 'root_bp_genomewide', 'n_shoot_hs_genomewide', 'V5', 'shoot_bp_genomewide')

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
          
largest5=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),5)))))
largest5
ind$largest5=ind$fam %in% names(largest5)                               


te=ind


get_largest_quantile_backgroundbox=function(feat){
 library(plyr)
 sub2=cbind(ind[ind$largest10 & !duplicated(ind$TEID),] %>% select_(feat) , ind[ind$largest10 & !duplicated(ind$TEID),] %>% select(fam))
 sub3=sub2 %>% group_by(fam) %>% summarize_all(funs(median(., na.rm=TRUE), min=quantile(., na.rm=TRUE, 0.25), max=quantile(.,na.rm=TRUE, 0.75)))
 sub3$sup=substr(sub3$fam,1,3)
 sub3$value=NA
 sub3$plt='point'
 sub1=cbind(ind %>% select_(feat) , ind %>% select(sup))
 af=sub1 %>% group_by(sup) %>% summarize_all(funs(median(., na.rm=TRUE), min=quantile(., na.rm=TRUE, 0.25), max=quantile(.,na.rm=TRUE, 0.75)))
 af$fam=af$sup
 af$value=NA
 af$plt='suppoint'
 colnames(af)[1]='fam'
 colnames(af)[5]='sup'
 colnames(af)[2:4]=paste0(colnames(af)[2:4], '_sup')
 sub3$median_sup=as.numeric(as.character(mapvalues(sub3$sup, from=af$sup, to=af$median_sup)))
 sub3$min_sup=as.numeric(as.character(mapvalues(sub3$sup, from=af$sup, to=af$max_sup)))
 sub3$max_sup=as.numeric(as.character(mapvalues(sub3$sup, from=af$sup, to=af$min_sup)))
 sub3=sub3[match(names(largest10), sub3$fam),]
 sub3$x=1:nrow(sub3)
 sub3$px=sub3$x
 sub3$x1=sub3$x+1
 sub3$px1=sub3$x1
 sub3$x[!duplicated(sub3$sup)]=sub3$x1[!duplicated(sub3$sup)]-0.75
 sub3$x1[which(!duplicated(sub3$sup))[-1]-1]=sub3$x1[which(!duplicated(sub3$sup))[-1]-1]-0.75

# sub4=merge(sub3, af, all=T)
 return(sub3)
}

## figure out percentages
get_largest_percents_backgroundbox=function(feat, invert=FALSE){
 library(plyr)
 sub2=cbind(ind[ind$largest10 & !duplicated(ind$TEID),] %>% select_(feat) , ind[ind$largest10 & !duplicated(ind$TEID),] %>% select(fam))
 classes=names(table(sub2[,feat, with=F]))
 if(invert){classes=rev(classes)}
 print(classes)
 if (length(classes)>2){
  print("warning - make sure the first class here is one you want to compare to all others")
  }
 sub3=sub2 %>% group_by(fam) %>% dplyr::summarize_all(funs(propFirst=sum(.==classes[1], na.rm=T)/sum(!is.na(.))))
 sub3$sup=substr(sub3$fam,1,3)
 sub3$plot='family'
 sub1=cbind(ind %>% select_(feat) , ind %>% select(sup))
 af=sub1 %>% group_by(sup) %>% dplyr::summarize_all(funs(propFirst=sum(.==classes[1], na.rm=T)/sum(!is.na(.))))
 af$fam=af$sup
 af$plot='superfamily'
 sub3$supperc=as.numeric(mapvalues(sub3$sup, from=af$sup, to=af$propFirst))
 sub3=sub3[match(names(largest10), sub3$fam),]
 sub3$x=1:nrow(sub3)
 sub3$px=sub3$x
 sub3$x1=sub3$x+1
 sub3$px1=sub3$x1
 sub3$x[!duplicated(sub3$sup)]=sub3$x1[!duplicated(sub3$sup)]-0.75
 sub3$x1[which(!duplicated(sub3$sup))[-1]-1]=sub3$x1[which(!duplicated(sub3$sup))[-1]-1]-0.75

# sub4=merge(sub3, af, all=T)
 return(sub3)
}


plot_percentages=function(feat, ylab='', invert=FALSE){
 ggplot(get_largest_percents_backgroundbox(feat, invert), aes(x=px, y=propFirst, fill=sup)) + 
                     geom_point(aes(color=sup), size=2) +
                      geom_col(aes(y=supperc), alpha=0.3, width=1) + 
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) +#ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
 }


### repeat for other figure of te vs flank


plotlargest=function(feat, ylab='', hline=-1000){
 ggplot(get_largest_quantile_backgroundbox(feat), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup)) + 
                     geom_hline(yintercept=hline, linetype='dashed', color='grey') +
                     geom_pointrange(fatten=4/3, size=1.5) + 
                     geom_rect(aes(xmin=x, xmax=x1, fill=sup, ymin=min_sup, ymax=max_sup), alpha=0.2, colour=NA) +
                     geom_point(aes(x=px+0.5, color=sup, y=median_sup), alpha=0.5, shape="-", size=1.5) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) + #ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
 }

pdf(paste0('figure5.', Sys.Date(), '.pdf'), 22,10)
gc=plotlargest('percGC', '% GC', hline=genomewide$percGC)
cg=plotlargest('nCG', 'Proportion\nCG methylatable', hline=genomewide$nCG)
mnase=plotlargest('shoot_prop', 'Proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen)
diversity=plotlargest('segsites.bp', 'Proportion\nsegregating sites', hline=genomewide.ss$genomewide_bp_segsites/genomewide$seqlen)
                              
gc.flank=plotlargest('percGC_1kbflank', 'Flanking % GC', hline=genomewide$percGC)
cg.flank=plotlargest('nCG_1kbflank', 'Flanking proportion\nCG methylatable', hline=genomewide$nCG)
mnase.flank=plotlargest('flank_shoot_prop', 'Flanking proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen)
diversity.flank=plotlargest('flank_segsites.bp', 'Flanking proportion\nsegregating sites', hline=genomewide.ss$genomewide_bp_segsites/genomewide$seqlen)

plot_grid(gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')

dev.off()
                              
                              
## supplemental base composition
pdf(paste0('pointrange_basecomp_flank.', Sys.Date(), '.pdf'), 22,14)
gc=plotlargest('percGC', '% GC', hline=genomewide$percGC)
cg=plotlargest('nCG', 'Proportion\nCG methylatable', hline=genomewide$nCG)
chg=plotlargest('nCHG', 'Proportion\nCHG methylatable', hline=genomewide$nCHG)
chh=plotlargest('nCHH', 'Proportion\nCHH methylatable', hline=genomewide$nCHH)
mnase.r=plotlargest('root_prop', 'Proportion Mnase\nhypersensitive (root)', hline=genomewide.mnase$root_bp_genomewide/genomewide$seqlen)
mnase.s=plotlargest('shoot_prop', 'Proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen)

gc.flank=plotlargest('percGC_1kbflank', 'Flanking % GC', hline=genomewide$percGC)
cg.flank=plotlargest('nCG_1kbflank', 'Flanking proportion\nCG methylatable', hline=genomewide$nCG)
chg.flank=plotlargest('nCHG_1kbflank', 'Flanking proportion\nCHG methylatable', hline=genomewide$nCHG)
chh.flank=plotlargest('nCHH_1kbflank', 'Flanking proportion\nCHH methylatable', hline=genomewide$nCHH)
mnase.flank.r=plotlargest('flank_root_prop', 'Flanking proportion Mnase\nhypersensitive (root)', hline=genomewide.mnase$root_bp_genomewide/genomewide$seqlen)
mnase.flank.s=plotlargest('flank_shoot_prop', 'Flanking proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen)

plot_grid(gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
          chg + ylim(0,0.12), chg.flank+ ylim(0,0.12), 
          chh+ ylim(0,0.2), chh.flank+ ylim(0,0.2) , 
          mnase.r+ ylim(0,0.15), mnase.flank.r + ylim(0,0.15),
          mnase.s + ylim(0,0.4), mnase.flank.s + ylim(0,0.4),
          labels = "AUTO", ncol = 2, align = 'v')
#plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(tel, cl, ingene, piece, disr + scale_x_discrete(labels=substr(names(largest10),1,3)[!duplicated(substr(names(largest10),1,3))]),  labels = "AUTO", ncol = 1, align = 'v')
## version with a legend.
#legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup))+ geom_pointrange(size=1)+ 
#                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
#plots <- plot_grid(tel, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plots <- plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(plots,legend, ncol = 2, align = 'v',  rel_widths = c(1, .1))                              
dev.off()
