library(cowplot)
library(data.table)
library(dplyr)
library(ggpubr)
library(scales)

source('../GenomeInfo.R')
source('color_palette.R')

tbl=fread('../te_age/tbl_age/B73_terminalbranchlength2018-10-27.txt')
ltr=fread('../te_age/B73v4_recovered_ages.txt')
ltr$TEID=gsub('B73v4', 'Zm00001d', ltr$tename)
ltr$tename=NULL

## always bring this one in!!!
allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
ind=merge(tbl, ltr, all=T)
ind$fam=substr(ind$TEID,1,8)
ind$sup=substr(ind$TEID,1,3)
ind=merge(ind, allte, all=T, by=c('TEID', 'fam', 'sup')) ## this leads to fam.x and fam.y, because the short ones we remove in the next line couldn't be built into trees!
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

## convert to real numbers
ind$age=ind$tbl
ind$age[!is.na(ind$k2p)]=ind$k2p[!is.na(ind$k2p)]
ind$mya=ind$age/3.3e-8/2/1e6
ind$ya=ind$age/3.3e-8/2
ind$tblmya=ind$tbl/3.3e-8/2/1e6

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
largest5=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),5)))))
largest5
ind$largest5=ind$fam %in% names(largest5)                               


te=ind                           
supAgeFun=function(superfam){
    ggplot(data.frame(te)[te$largest10 & te$sup==superfam & !is.na(te$sup),], aes(x=mya*1e6, fill=sup, group=fam)) + geom_histogram(binwidth=1e4) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Age (million years)', breaks=c(0,1e6,2e6), labels=c(0,1,2), limits=c(0,2.1e6))

}
supAgeFun5=function(superfam){
    ggplot(data.frame(te)[ te$largest5 & te$sup==superfam & !is.na(te$sup),], aes(x=mya*1e6, fill=sup, group=fam)) + geom_histogram(binwidth=1e4) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10)) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Age (million years)', breaks=c(0,1e6,2e6), labels=c(0,1,2), limits=c(0,2.1e6)) +
           scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}
supAgeFuntbl=function(superfam){
    ggplot(data.frame(te)[te$largest10 & te$sup==superfam & !is.na(te$sup),], aes(x=tblmya*1e6, fill=sup, group=fam)) + geom_histogram(binwidth=1e4) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Age (million years)', breaks=c(0,1e6,2e6), labels=c(0,1,2), limits=c(0,2.1e6))

}
                              
supAgeFun1my=function(superfam){
    ggplot(data.frame(te)[te$largest10 & te$sup==superfam & !is.na(te$sup),], aes(x=mya*1e6, fill=sup, group=fam)) + geom_histogram(binwidth=1e4) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Age (million years)', breaks=c(0,5e5,1e6), labels=c(0,0.5,1), limits=c(0,1.01e6))

}
supAgeFun51my=function(superfam){
    ggplot(data.frame(te)[ te$largest5 & te$sup==superfam & !is.na(te$sup),], aes(x=mya*1e6, fill=sup, group=fam)) + geom_histogram(binwidth=1e4) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10)) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Age (million years)', breaks=c(0,5e5,1e6), labels=c(0,0.5,1), limits=c(0,1.01e6)) +
           scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}
supAgeFuntbl1my=function(superfam){
    ggplot(data.frame(te)[te$largest10 & te$sup==superfam & !is.na(te$sup),], aes(x=tblmya*1e6, fill=sup, group=fam)) + geom_histogram(binwidth=1e4) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Age (million years)', breaks=c(0,5e5,1e6), labels=c(0,0.5,1), limits=c(0,1.01e6))

}
                              
                              

## figure 3
pdf(paste0('figure3.', Sys.Date(), '.pdf'), 10,10)
sups=ggplot(te[!is.na(te$sup),], aes(x=mya*1e6, fill=factor(sup, levels=TESUPFACTORLEVELS))) + geom_histogram(binwidth=1e4) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col, name='') + scale_x_continuous(name='Age (million years)', breaks=c(0,1e6,2e6), labels=c(0,1,2), limits=c(0,2.1e6)) + ylab('Number copies')
#DHH= ggplot(te[te$chr==1 & te$largest10 & te$famsize>10 & te$sup=='DHH',], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
#           theme(strip.background = element_blank(),strip.text.x = element_blank()) + scale_fill_manual(values=dd.col)
DHH = supAgeFun('DHH')
DTA = supAgeFun('DTA')
DTC = supAgeFun('DTC')
DTH = supAgeFun('DTH')
DTM = supAgeFun('DTM')                              
DTT = supAgeFun('DTT')                              
DTX = supAgeFun('DTX')
RLC = supAgeFun('RLC')                             
RLG = supAgeFun('RLG')                             
RLX = supAgeFun('RLX')                                                          
RIT = supAgeFun('RIT')                              
RIL = supAgeFun('RIL')                              
RST = supAgeFun('RST')                              

DHH5 = supAgeFun5('DHH')
DTA5 = supAgeFun5('DTA')
DTC5 = supAgeFun5('DTC')
DTH5 = supAgeFun5('DTH')
DTM5 = supAgeFun5('DTM')                              
DTT5 = supAgeFun5('DTT')                              
DTX5 = supAgeFun5('DTX')
RLC5 = supAgeFun5('RLC')                             
RLG5 = supAgeFun5('RLG')                             
RLX5 = supAgeFun5('RLX')                                                          
RIT5 = supAgeFun5('RIT')
RIL5 = supAgeFun5('RIL')
RST5 = supAgeFun5('RST')                              
                                     
                              
rightside=plot_grid(DHH, DTA, DTC, DTH, DTM, DTT, RLC, RLG, RLX, RST, labels=c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'), ncol=1)
bigones=plot_grid(DHH + theme(legend.position="none"), RLC + theme(legend.position="none"), RLG + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')
almostbig=plot_grid(DHH5+ theme(legend.position="none"), DTA5+ theme(legend.position="none"), RLC5+ theme(legend.position="none"), RLG5+ theme(legend.position="none"), labels=c('B', 'C', 'D', 'E', 'F', 'G'), ncol=1)
bigones5=plot_grid(DHH5 + theme(legend.position="none"), RLC5 + theme(legend.position="none"), RLG5 + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')

#plot_grid(sups, rightside, labels=c('A', ''), ncol=2, align='v')
plot_grid(sups, bigones, labels=c('A', ''), ncol=2, align='v')

plot_grid(sups, almostbig, labels=c('A', ''), ncol=2, align='v')

plot_grid(sups, bigones5, labels=c('A', ''), ncol=2, align='h', rel_widths = c(1.5, 1.1))

## also try 1 mya limit
sups1my=ggplot(te[!is.na(te$sup),], aes(x=mya*1e6, fill=factor(sup, levels=TESUPFACTORLEVELS))) + geom_histogram(binwidth=1e4) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col, name='') + scale_x_continuous(name='Age (million years)', breaks=c(0,5e5,1e6), labels=c(0,0.5,1), limits=c(0,1.01e6)) + ylab('Number copies')

#rightside1my=plot_grid(DHH + xlim(0,1e6), DTA + xlim(0,1e6), DTC + xlim(0,1e6), DTH + xlim(0,1e6), DTM + xlim(0,1e6), DTT + xlim(0,1e6), RLC + xlim(0,1e6), RLG + xlim(0,1e6), RLX + xlim(0,1e6), RST + xlim(0,1e6), labels=c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'), ncol=1)
#bigones1my=plot_grid(DHH + theme(legend.position="none"), RLC + theme(legend.position="none"), RLG + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')

DHH51my=supAgeFun51my('DHH')           
DTA51my=supAgeFun51my('DTA')
RLC51my=supAgeFun51my('RLC')
RLG51my=supAgeFun51my('RLG')
almostbig1my=plot_grid(DHH51my+ theme(legend.position="none"), DTA51my + theme(legend.position="none"), RLC51my+ theme(legend.position="none"), RLG51my+ theme(legend.position="none"), labels=c('B', 'C', 'D', 'E'), ncol=1)
#bigones51my=plot_grid(DHH5 + theme(legend.position="none"), RLC5 + theme(legend.position="none"), RLG5 + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')

#plot_grid(sups, rightside, labels=c('A', ''), ncol=2, align='v')
#plot_grid(sups + xlim(0,1e6), bigones1my, labels=c('A', ''), ncol=2, align='v')

plot_grid(sups1my, almostbig1my, labels=c('A', ''), ncol=2, align='v')

#plot_grid(sups+ xlim(0,1e6), bigones51my, labels=c('A', ''), ncol=2, align='h', rel_widths = c(1.5, 1.1))
                              
                              
  
dev.off()
                              
pdf('supplemental_TEage.pdf', 10,15)
## all sups, 10 largest fams
plot_grid(DHH+ theme(legend.position="none", axis.title.x=element_blank()), DTA+ theme(legend.position="none", axis.title.x=element_blank()), DTC+ theme(legend.position="none", axis.title.x=element_blank()),
          DTH+ theme(legend.position="none", axis.title.x=element_blank()), DTM+ theme(legend.position="none", axis.title.x=element_blank()), DTT+ theme(legend.position="none", axis.title.x=element_blank()),
          DTX+ theme(legend.position="none", axis.title.x=element_blank()), RLC+ theme(legend.position="none", axis.title.x=element_blank()), RLG+ theme(legend.position="none", axis.title.x=element_blank()),
          RLX+ theme(legend.position="none", axis.title.x=element_blank()), RIT+ theme(legend.position="none"), RIL+ theme(legend.position="none"), 
          RST+ theme(legend.position="none"), labels="AUTO", ncol=3)
## switch main plots to tbl age for LTR's
RLCtbl = supAgeFuntbl('RLC')                             
RLGtbl = supAgeFuntbl('RLG')                             
RLXtbl = supAgeFuntbl('RLX')                                                          

ltrssup=ggplot(data.frame(te)[!is.na(te$sup) & te$sup %in% c('RLC', 'RLG', 'RLX'),], aes(x=mya*1e6, fill=factor(sup, levels=TESUPFACTORLEVELS))) + geom_histogram(binwidth=1e4) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col, name='') + scale_x_continuous(name='Age (million years)', breaks=c(0,1e6,2e6), labels=c(0,1,2), limits=c(0,2.1e6)) + ylab('Number copies')
ltrssuptbl=ggplot(data.frame(te)[!is.na(te$sup) & te$sup %in% c('RLC', 'RLG', 'RLX'),], aes(x=tblmya*1e6, fill=factor(sup, levels=TESUPFACTORLEVELS))) + geom_histogram(binwidth=1e4) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col, name='') + scale_x_continuous(name='Age (million years)', breaks=c(0,1e6,2e6), labels=c(0,1,2), limits=c(0,2.1e6)) + ylab('Number copies')
ltrs=plot_grid(ltrssup, RLC + theme(legend.position="none"), RLG + theme(legend.position="none"), RLX + theme(legend.position="none"), labels="AUTO", ncol=1, align='h')
ltrstbl=plot_grid(ltrssuptbl, RLCtbl + theme(legend.position="none"), RLGtbl + theme(legend.position="none"), RLXtbl + theme(legend.position="none"), labels=c('E', 'F', 'G', 'H'), ncol=1, align='h')
plot_grid(ltrs, ltrstbl, labels='', ncol=2, align='v')
dev.off()
                              
pdf('supplemental_LTR_TBL.pdf', 15,5)
## compare LTR vs tbl age! - spearman corr coef on plot
ggplot(data.frame(te)[te$sup %in% c('RLC', 'RLG', 'RLX'),], aes(x=tblmya, y=mya, color=sup)) + geom_point(alpha=0.1) + 
                              scale_color_manual(name='', values=dd.col) + xlab('Terminal branch length age (million years)') + 
                              ylab('LTR-LTR age (million years)') + xlim(0,2.1) + ylim(0,2.1) + facet_wrap(~sup, ncol=3) +
                              stat_cor(aes(color = sup), label.x = 1, method = "spearman")

dev.off()                          
                              
                              
