library(RColorBrewer)
library(treemapify)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(plyr)
library(dplyr)
library(scales)
library(grid)

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
te$famFactor=factor(te$fam, levels=names(largest5))

supChrFun=function(superfam){
    ggplot(te[te$chr==1 & te$sup==superfam,], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), 
#                 strip.text.y = element_text(angle = 360), 
                 strip.text.y=element_blank(),
                 axis.text=element_text(size=10)) + 
          scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) +
          scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}
supChrFun5=function(superfam, chrom=1){
    ggplot(te[te$chr==chrom & te$largest5 & te$sup==superfam,], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), 
#                 strip.text.y = element_text(angle = 360), 
                 strip.text.y=element_blank(),
                 axis.text=element_text(size=10)) + 
          scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) +
          scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}
                              
supChrFun5Label=function(superfam, chrom=1){
    ggplot(te[te$chr==chrom & te$largest5 & te$sup==superfam,], aes(x=start, fill=sup, group=famFactor)) + geom_histogram(binwidth=1e6) + facet_wrap(~famFactor, ncol=1, scales='free_y', drop=T, strip.position='right')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), 
                 strip.text.y = element_text(angle = 360), 
 #                strip.text.y=element_blank(),
                 axis.text=element_text(size=10)) + 
          scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) +
          scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}
                              
supChrFun5multsups=function(superfam){
    ggplot(te[te$chr==1 & te$largest5 & te$sup%in%superfam,], aes(x=start, fill=sup, group=fam)) + geom_histogram(binwidth=1e6) + facet_wrap(~fam, ncol=1, drop=T, strip.position='right', scales='free_y')+ 
           theme(strip.background = element_blank(),strip.text.x = element_blank(), 
#                 strip.text.y = element_text(angle = 360), 
                 strip.text.y=element_blank(),
                 axis.text=element_text(size=10)) + 
          scale_fill_manual(values=dd.col) +
           ylab('') + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + 
          scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))

}                              

pdf(paste0('figure3.chromosome1_dists.', Sys.Date(), '.pdf'), 10,10)
sups=ggplot(te[te$chr==1,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10), legend.position='bottom') +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('Number copies') +  scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))
supsNamed=ggplot(te[te$chr==1,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + 
                              facet_wrap(factor(sup, levels=TESUPFACTORLEVELS)~., ncol=1, scales='free_y', strip.position='right')+ 
                              theme(strip.background = element_blank(),strip.text = element_text(angle = 270), axis.text=element_text(size=10), legend.position='none') +  
                              scale_fill_manual(values=dd.col) + 
                              scale_x_continuous(name='Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + 
                              ylab('Number copies') +  scale_y_continuous(breaks=scales::pretty_breaks(2), limits=c(0,NA))
                              
                              
## have to go into the weeds to get each facet text to be a different color!
g <- ggplot_gtable(ggplot_build(supsNamed))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- dd.col
k <- 1
for (i in stripr) {
#  j <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  sup=as.character(g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label)
  print(sup)
  g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- fills[sup]
#  k <- k+1
}
#grid.draw(g)
                              
                              
                              
                 
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
              
#DHHDTTRLCRLG5=supChrFun5multsups(c('DHH', 'DTT', 'RLC', 'RLG'))
                              
#rightside=plot_grid(DHH, DTA, DTC, DTH, DTM, DTT, DTX, RLC, RLG, RLX, RIL, RIT, RST, labels=c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N'), ncol=1)
#bigones=plot_grid(DHH + theme(legend.position="none"), RLC + theme(legend.position="none"), RLG + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')
almostbigH=plot_grid(DHH5+ theme(legend.position="none", axis.text.x=element_blank(), axis.title.x=element_blank()), DTH5+ theme(legend.position="none", axis.text.x=element_blank(), axis.title.x=element_blank()), RLC5+ theme(legend.position="none", axis.text.x=element_blank(), axis.title.x=element_blank()), RLG5+ theme(legend.position="none"), labels=c('B', 'C', 'D', 'E', 'F', 'G'), ncol=1)
almostbigT=plot_grid(DHH5+ theme(legend.position="none", axis.text.x=element_blank(), 
                                 axis.title.x=element_blank()), 
                     DTT5+ theme(legend.position="none", axis.text.x=element_blank(), 
                                 axis.title.x=element_blank()), 
                     RLC5+ theme(legend.position="none", axis.text.x=element_blank(), 
                                 axis.title.x=element_blank()), 
                     RLG5+ theme(legend.position="none"), labels=c('B', 'C', 'D', 'E', 'F', 'G'), ncol=1)
bigones5=plot_grid(DHH5 + theme(legend.position="none"), RLC5 + theme(legend.position="none"), RLG5 + theme(legend.position="none"), labels=c('B', 'C', 'D'), ncol=1, align='h')
foursups=plot_grid(DHHDTTRLCRLG5 + theme(legend.position='none'), labels=c('B'), ncol=1, align='h')
#plot_grid(sups, rightside, labels=c('A', ''), ncol=2, align='v')
#plot_grid(sups+ theme(legend.position="bottom"), bigones, labels=c('A', ''), ncol=2, align='v')

plot_grid(sups+ theme(legend.position="bottom"), almostbigH, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))
plot_grid(sups+ theme(legend.position="bottom"), almostbigT, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))

plot_grid(sups+ theme(legend.position="bottom"), foursups, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))
             
plot_grid(sups+ theme(legend.position="bottom"), bigones5, labels=c('A', ''), ncol=2, align='h', rel_widths = c(1.5, 1.1))

                              
plot_grid(g, almostbigT, labels=c('A', ''), ncol=2, align='v', rel_widths = c(1.5, 1.1))
plot_grid(g, bigones5, labels=c('A', ''), ncol=2, align='h', rel_widths = c(1.5, 1.1))

                              
dev.off()

pdf(paste0('all_chroms_supp.', Sys.Date(), '.pdf'), 20,10)
sups1=ggplot(te[te$chr==1,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 1 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('Number copies') + theme(legend.position="none")
sups2=ggplot(te[te$chr==2,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 2 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups3=ggplot(te[te$chr==3,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 3 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups4=ggplot(te[te$chr==4,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 4 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups5=ggplot(te[te$chr==5,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 5 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups6=ggplot(te[te$chr==6,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 6 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups7=ggplot(te[te$chr==7,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 7 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups8=ggplot(te[te$chr==8,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 8 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups9=ggplot(te[te$chr==9,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 9 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')+ theme(legend.position="none")
sups10=ggplot(te[te$chr==10,], aes(x=start, fill=sup)) + geom_histogram(binwidth=1e6) + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), ncol=1, scales='free_y')+ theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text=element_text(size=10)) +  scale_fill_manual(values=dd.col) + scale_x_continuous(name='Chr 10 Position (Mb)', breaks=c(0,1e8,2e8, 3e8), labels=c(0,100,200,300)) + ylab('')

                              
plot_grid(sups1, sups2, sups3, sups4, sups5, sups6, sups7, sups8, sups9, sups10, ncol=10, align='h', rel_widths=c(1,1,1,1,1,1,1,1,1,1.5))
dev.off()
                              
                              
### to tif format, this is S2_Fig.tif
                              
tiff(paste0('all_chroms_supp.', Sys.Date(), '.tif'), 20,10, res=300, units='in')
plot_grid(sups1, sups2, sups3, sups4, sups5, sups6, sups7, sups8, sups9, sups10, ncol=10, align='h', rel_widths=c(1,1,1,1,1,1,1,1,1,1.5))

dev.off()
                              
pdf(paste0('supp_chromsFam.', Sys.Date(), '.pdf'), 10,10)
plot_grid(DHH5+ theme(legend.position="none"), DTA5+ theme(legend.position="none"),DTC5+ theme(legend.position="none"),DTH5+ theme(legend.position="none"),DTM5+ theme(legend.position="none"),DTT5+ theme(legend.position="none"),DTX5+ theme(legend.position="none"),RLC5+ theme(legend.position="none"),RLG5+ theme(legend.position="none"),RLX5+ theme(legend.position="none"),RIL5+ theme(legend.position="none"),RIT5+ theme(legend.position="none"),RST5+ theme(legend.position="none"), labels='AUTO', ncol=4, align='h')

                              
## weird, some DNA tes not on a chromosome?!?!?
                              
for(chrom in 2:10){
          DHH5 = supChrFun5('DHH', chrom=chrom)
          DTA5 = supChrFun5('DTA', chrom=chrom)
          DTC5 = supChrFun5('DTC', chrom=chrom)
          DTH5 = supChrFun5('DTH', chrom=chrom)
          DTM5 = supChrFun5('DTM', chrom=chrom)                              
          DTT5 = supChrFun5('DTT', chrom=chrom)                              
          DTX5 = supChrFun5('DTX', chrom=chrom)
          RLC5 = supChrFun5('RLC', chrom=chrom)                             
          RLG5 = supChrFun5('RLG', chrom=chrom)                             
          RLX5 = supChrFun5('RLX', chrom=chrom)                                                          
          RST5 = supChrFun5('RST', chrom=chrom)                              
          RIL5 = supChrFun5('RIL', chrom=chrom)                              
          RIT5 = supChrFun5('RIT', chrom=chrom)                              
          print(plot_grid(DHH5+ theme(legend.position="none"), DTA5+ theme(legend.position="none"),DTC5+ theme(legend.position="none"),DTH5+ theme(legend.position="none"),DTM5+ theme(legend.position="none"),DTT5+ theme(legend.position="none"),DTX5+ theme(legend.position="none"),RLC5+ theme(legend.position="none"),RLG5+ theme(legend.position="none"),RLX5+ theme(legend.position="none"),RIL5+ theme(legend.position="none"),RIT5+ theme(legend.position="none"),RST5+ theme(legend.position="none"), labels='AUTO', ncol=4, align='h'))
}         
dev.off()                              

### to tif format, this is S3_Fig.tif
                              
DHH5 = supChrFun5Label('DHH')
DTA5 = supChrFun5Label('DTA')
DTC5 = supChrFun5Label('DTC')
DTH5 = supChrFun5Label('DTH')
DTM5 = supChrFun5Label('DTM')                              
DTT5 = supChrFun5Label('DTT')                              
DTX5 = supChrFun5Label('DTX')
RLC5 = supChrFun5Label('RLC')                             
RLG5 = supChrFun5Label('RLG')                             
RLX5 = supChrFun5Label('RLX')                                                          
RST5 = supChrFun5Label('RST')                              
RIL5 = supChrFun5Label('RIL')                              
RIT5 = supChrFun5Label('RIT') 
                              
tiff(paste0('supp_chromsFam.', Sys.Date(), '.tif'), 14,10, res=300, units='in')
plot_grid(DHH5+ theme(legend.position="none"), DTA5+ theme(legend.position="none"),DTC5+ theme(legend.position="none"),DTH5+ theme(legend.position="none"),DTM5+ theme(legend.position="none"),DTT5+ theme(legend.position="none"),DTX5+ theme(legend.position="none"),RLC5+ theme(legend.position="none"),RLG5+ theme(legend.position="none"),RLX5+ theme(legend.position="none"),RIL5+ theme(legend.position="none"),RIT5+ theme(legend.position="none"),RST5+ theme(legend.position="none"), labels='AUTO', ncol=4, align='h')

dev.off()
                              
                              
## see chr dist here!
sapply(names(largest5), function(x) table(ind$chr[ind$fam==x]))
                              
                              
                              
                              
