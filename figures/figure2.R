library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)
library(scales)

source('../GenomeInfo.R')
source('color_palette.R')

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
colnames(gene)[3]='TEID'

ind=merge(techar, gene, all=T)
ind$ingene=ind$closest==0
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

nrow(ind)
nrow(techar)
nrow(gene)

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
          
largest5=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),5)))))
largest5
ind$largest5=ind$fam %in% names(largest5)                               


te=ind

supChrFun=function(superfam){
    ggplot(te[te$chr==1 & te$sup==superfam,], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10)) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300))

}
supChrFun5=function(superfam){
    ggplot(te[te$chr==1 & te$largest5 & te$sup==superfam,], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10)) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) +
          scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}
                              
supChrFun5multsups=function(superfam){
    ggplot(te[te$chr==1 & te$largest5 & te$sup%in%superfam,], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, drop=T, strip.position='right', scales='free_y')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.text=element_text(size=10)) + scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + 
          scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}                              

pdf(paste0('figure2.chromosome1_dists.', Sys.Date(), '.pdf'), 10,10)
sups=ggplot(te[te$chr==1,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('Number copies')
#DHH= ggplot(te[te$chr==1 & te$largest10 & te$famsize>10 & te$sup=='DHH',], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
#           theme(strip.background = element_blank(),strip.text.x = element_blank()) + scale_fill_manual(values=dd.col)
DHH = supChrFun('DHH')
DTA = supChrFun('DTA')
DTC = supChrFun('DTC')
DTH = supChrFun('DTH')
DTM = supChrFun('DTM')                              
DTT = supChrFun('DTT')                              
DTX = supChrFun('DTX')
RLC = supChrFun('RLC')                             
RLG = supChrFun('RLG')                             
RLX = supChrFun('RLX')                                                          
RST = supChrFun('RST')
RIL = supChrFun('RIL')
RIT = supChrFun('RIT')

DHH5 = supChrFun5('DHH')
DTA5 = supChrFun5('DTA')
DTC5 = supChrFun5('DTC')
DTH5 = supChrFun5('DTH')
DTM5 = supChrFun5('DTM')                              
DTT5 = supChrFun5('DTT')                              
DTX5 = supChrFun5('DTX')
RLC5 = supChrFun5('RLC')                             
RLG5 = supChrFun5('RLG')                             
RLX5 = supChrFun5('RLX')                                                          
RST5 = supChrFun5('RST')                              
RIL5 = supChrFun5('RIL')                              
RIT5 = supChrFun5('RIT')                              
              
DHHDTTRLCRLG5=supChrFun5multsups(c('DHH', 'DTT', 'RLC', 'RLG'))
                              
rightside=plot_grid(DHH, DTA, DTC, DTH, DTM, DTT, DTX, RLC, RLG, RLX, RIL, RIT, RST, labels=c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N'), ncol=1)
bigones=plot_grid(DHH + theme(legend.position="none"), RLC + theme(legend.position="none"), RLG + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')
almostbigH=plot_grid(DHH5+ theme(legend.position="none"), DTH5+ theme(legend.position="none"), RLC5+ theme(legend.position="none"), RLG5+ theme(legend.position="none"), labels=c('B', 'C', 'D', 'E', 'F', 'G'), ncol=1)
almostbigT=plot_grid(DHH5+ theme(legend.position="none"), DTT5+ theme(legend.position="none"), RLC5+ theme(legend.position="none"), RLG5+ theme(legend.position="none"), labels=c('B', 'C', 'D', 'E', 'F', 'G'), ncol=1)
bigones5=plot_grid(DHH5 + theme(legend.position="none"), RLC5 + theme(legend.position="none"), RLG5 + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')
foursups=plot_grid(DHHDTTRLCRLG5 + theme(legend.position='none'), labels=c('B'), ncol=1, align='h')
#plot_grid(sups, rightside, labels=c('A', ''), ncol=2, align='v')
plot_grid(sups, bigones, labels=c('A', ''), ncol=2, align='v')

plot_grid(sups, almostbigH, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))
plot_grid(sups, almostbigT, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))

plot_grid(sups, foursups, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))
             
plot_grid(sups, bigones5, labels=c('A', ''), ncol=2, align='h', rel_widths = c(1.5, 1.1))

dev.off()
