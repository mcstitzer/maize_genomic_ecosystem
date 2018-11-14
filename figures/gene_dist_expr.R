library(rtracklayer)
library(stringr)
library(lme4)
library(dplyr)
library(reshape2)
library(cowplot)
library(data.table)

source('../GenomeInfo.R')
source('color_palette.R')

tgem=fread('../genes/B73_closest_gene_expression.maizegdbWalley.txt')

colnames(tgem)[3]='TEID'
allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
ind=merge(tgem, allte, all=T)
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
ind$famsize=as.numeric(table(ind$fam)[ind$fam])
ind=ind[,-c('chr', 'start', 'end', 'strand', 'source', 'type', 'tebp', 'tespan', 'pieces', 'disruptor', 'largest10', 'geneID')]

### try out normalizing by median, max, and transforming to stdnormal                            
ind$gene_max=apply(ind[,7:29], 1, function(x) max(x, na.rm=T))
ind$gene_max[is.infinite(ind$gene_max)]=NA ## because if all NA, this becomes -Inf
                   
ind.scaledmax=ind
ind.scaledmax[,7:29]=ind[,7:29]/ind$gene_max

ind$gene_mean=apply(ind[,7:29], 1, function(x) mean(x, na.rm=T))
ind$gene_mean[is.infinite(ind$gene_mean)]=NA ## because if all NA, this becomes -Inf
ind$gene_sd=apply(ind[,7:29], 1, function(x) sd(x, na.rm=T))
ind$gene_sd[is.infinite(ind$gene_sd)]=NA ## because if all NA, this becomes -Inf

ind.standardized=ind
ind.standardized[,7:29]=(ind[,7:29]-ind$gene_mean)/ind$gene_sd
                  
                  
ge=melt(ind, id.vars=c("sup", 'fam', 'TEID', 'closest', 'closestgene', 'closestgenetype', 'famsize'))
ge$bins=cut(ge$closest, breaks=seq(0,750000, length.out=51), include.lowest=T, labels=seq(0,750000, length.out=51)[-1])
ge$bins.10kb=cut(ge$closest, breaks=seq(0,10000, length.out=21), include.lowest=T, labels=seq(0,10000, length.out=21)[-1])
ge$bins.20kb=cut(ge$closest, breaks=seq(0,20000, length.out=41), include.lowest=T, labels=seq(0,20000, length.out=41)[-1])
ge$bins.10kb100bp=cut(ge$closest, breaks=seq(0,10000, length.out=101), include.lowest=T, labels=seq(0,10000, length.out=101)[-1])

#supmed=tgem %>% group_by(sup) %>% summarize_if(.predicate=function(x) is.numeric(x), .funs=funs(median(., na.rm=T)))
#supge=melt(supmed, id.vars=c("sup", 'closest','famsize'))

fammed=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize)  %>% dplyr::summarize(value=median(value, na.rm=T))
fammean=ge %>% group_by(sup, fam, bins, bins.10kb, bins.20kb, bins.10kb100bp, variable) %>% dplyr::summarize(value=mean(value, na.rm=T))
#famge=melt(fammed, id.vars=c("sup", 'fam', 'closest','famsize'))

#d.lg=d.l %>% group_by(sup, fam, variable, distance, context) %>% dplyr::summarize(value=mean(value, na.rm=T))
#d.mg=d.m%>% group_by(sup, fam, variable, distance, context, famsize) %>% dplyr::summarize(value=mean(value, na.rm=T))
  

#pdf('gene_dist_expr.pdf', 20,20)
#ggplot(supge, aes(x=log10(closest), y=log10(value), color=sup)) + geom_point()  + scale_color_manual(values=dd.col)
#ggplot(famge, aes(x=log10(closest), y=log10(value), color=factor(variable))) + geom_point()  + scale_color_manual(values=dd.col)
#ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col)
#ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#    
#ggplot(ge, aes(x=bins.10kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge, aes(x=bins.10kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge, aes(x=bins.10kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#
#                                                    
#ggplot(ge, aes(x=bins.20kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge, aes(x=bins.20kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge, aes(x=bins.20kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#      
#ggplot(ge, aes(x=bins.10kb100bp, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge, aes(x=bins.10kb100bp, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge, aes(x=bins.10kb100bp, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#                                                          
#ggplot(ge[ge$famsize>=10,], aes(x=bins.10kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge[ge$famsize>=10,], aes(x=bins.10kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge[ge$famsize>=10,], aes(x=bins.10kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#
#ggplot(ge[ge$famsize>=10,], aes(x=bins.20kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge[ge$famsize>=10,], aes(x=bins.20kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge[ge$famsize>=10,], aes(x=bins.20kb, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#      
#ggplot(ge[ge$famsize>=10,], aes(x=bins.10kb100bp, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
#ggplot(ge[ge$famsize>=10,], aes(x=bins.10kb100bp, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
#ggplot(ge[ge$famsize>=10,], aes(x=bins.10kb100bp, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
#                                                       
##ggplot(ge, aes(x=bins, y=value, color=sup, group=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + stat_smooth(geom='loess', se=F, alpha=0.3)
##ggplot(ge, aes(x=bins, y=value, color=sup, group=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + stat_smooth(geom='loess', se=F, alpha=0.3)
#
#dev.off()
                                                                                               
                                                                                               
fammed$distance=as.numeric(as.character(fammed$bins.10kb100bp))                                                                                               
pdf('expression_decay.pdf', 25,25)
#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
tissue='anther'
## ten largest families of each sup
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank()) + ylim(0,100)
## all families >10
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,100))
#print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
#                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
#                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10, alpha by famsize
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,100))

                               
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=log10(value), color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank())+ ylim(0,2)
## all families >10
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=log10(value), col=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,2))
#print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
#                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
#                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10, alpha by famsize
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=log10(value), col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,2))
                                                                                            
dev.off()                                           

pdf('expression_decay.2kb.pdf', 25,25)
#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
tissue='anther'
## ten largest families of each sup
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank()) + ylim(0,100)+ xlim(0,2100)
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(0,100)+ xlim(0,2100))

ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=log10(value), color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(0,2)+ xlim(0,2100)
## all families >10, alpha by famsize
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=log10(value), col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(0,2)+ xlim(0,2100))
                                                                                            
dev.off()                                           


                  
                  
                  
                  
                  
                  
#### scaled max
ge=melt(ind.scaledmax, id.vars=c("sup", 'fam', 'TEID', 'closest', 'closestgene', 'closestgenetype', 'famsize'))
ge$bins=cut(ge$closest, breaks=seq(0,750000, length.out=51), include.lowest=T, labels=seq(0,750000, length.out=51)[-1])
ge$bins.10kb=cut(ge$closest, breaks=seq(0,10000, length.out=21), include.lowest=T, labels=seq(0,10000, length.out=21)[-1])
ge$bins.20kb=cut(ge$closest, breaks=seq(0,20000, length.out=41), include.lowest=T, labels=seq(0,20000, length.out=41)[-1])
ge$bins.10kb100bp=cut(ge$closest, breaks=seq(0,10000, length.out=101), include.lowest=T, labels=seq(0,10000, length.out=101)[-1])

fammed=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize)  %>% dplyr::summarize(value=median(value, na.rm=T))
fammean=ge %>% group_by(sup, fam, bins, bins.10kb, bins.20kb, bins.10kb100bp, variable) %>% dplyr::summarize(value=mean(value, na.rm=T))
## keep only those with at least one observation of a gene in each bin.
fammedEach=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize) %>% filter(all(seq(0,10000, length.out=101)[-1] %in% bins.10kb100bp)) %>% dplyr::summarize(value=median(value, na.rm=T))
## keep only those with at least one observation of a gene in each bin.
fammed10=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize) %>% filter(n() >= 10) %>% dplyr::summarize(value=median(value, na.rm=T))
               
fammedEach$distance=as.numeric(as.character(fammedEach$bins.10kb100bp))                   
fammed10$distance=as.numeric(as.character(fammed10$bins.10kb100bp))                   
                                                                                            
fammed$distance=as.numeric(as.character(fammed$bins.10kb100bp))                                                                                               
pdf('expression_decay.scaledmax.pdf', 25,25)
#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
tissue='anther'
## ten largest families of each sup
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank()) + ylim(0,1)
## all families >10
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,1))
#print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
#                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
#                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10, alpha by famsize
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,1))

                               
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=log10(value), color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank())+ ylim(0,2)
## all families >10
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=log10(value), col=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,2))
#print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
#                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
#                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10, alpha by famsize
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=log10(value), col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(0,2))
                                                                                            
dev.off()                                           

pdf('expression_decay.2kb.scaledmax.pdf', 25,25)
#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
tissue='anther'
## ten largest families of each sup
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank()) + ylim(0,1)+ xlim(0,2100)
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(0,1)+ xlim(0,2100))

ggplot(data.frame(fammed10[fammed10$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank()) + ylim(0,1)+ xlim(0,2100)
print(ggplot(fammed10[fammed10$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(0,1)+ xlim(0,2100))
                         
dev.off()                    
                 
                  
                  
                  
#### standardized 
ge=melt(ind.standardized, id.vars=c("sup", 'fam', 'TEID', 'closest', 'closestgene', 'closestgenetype', 'famsize'))
ge$bins=cut(ge$closest, breaks=seq(0,750000, length.out=51), include.lowest=T, labels=seq(0,750000, length.out=51)[-1])
ge$bins.10kb=cut(ge$closest, breaks=seq(0,10000, length.out=21), include.lowest=T, labels=seq(0,10000, length.out=21)[-1])
ge$bins.20kb=cut(ge$closest, breaks=seq(0,20000, length.out=41), include.lowest=T, labels=seq(0,20000, length.out=41)[-1])
ge$bins.10kb100bp=cut(ge$closest, breaks=seq(0,10000, length.out=101), include.lowest=T, labels=seq(0,10000, length.out=101)[-1])

fammed=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize)  %>% dplyr::summarize(value=median(value, na.rm=T))
fammean=ge %>% group_by(sup, fam, bins, bins.10kb, bins.20kb, bins.10kb100bp, variable) %>% dplyr::summarize(value=mean(value, na.rm=T))
## keep only those with at least one observation of a gene in each bin.
fammedEach=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize) %>% filter(all(seq(0,10000, length.out=101)[-1] %in% bins.10kb100bp)) %>% dplyr::summarize(value=median(value, na.rm=T))
## keep only those with at least one observation of a gene in each bin.
fammed10=ge %>% group_by(sup, fam, bins.10kb100bp, variable, famsize) %>% filter(n() >= 10) %>% dplyr::summarize(value=median(value, na.rm=T))
               
fammedEach$distance=as.numeric(as.character(fammedEach$bins.10kb100bp))                   
fammed10$distance=as.numeric(as.character(fammed10$bins.10kb100bp))                   
fammed$distance=as.numeric(as.character(fammed$bins.10kb100bp))                                                                                               
pdf('expression_decay.standardized.pdf', 25,25)
#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
tissue='anther'
## ten largest families of each sup
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank()) + ylim(-5,5)
## all families >10
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(-5,5))
#print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
#                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
#                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10, alpha by famsize
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.x=element_blank())+ ylim(-5,5))

                               
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=log10(value), color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(variable~sup, ncol=13) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank())+ ylim(0,2)
                                                                                           
dev.off()                                           

pdf('expression_decay.2kb.standardized.pdf', 25,25)
#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
## ten largest families of each sup
ggplot(data.frame(fammed[fammed$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank()) + ylim(-5,5)+ xlim(0,2100)
print(ggplot(fammed[fammed$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(-5,5)+ xlim(0,2100))
ggplot(data.frame(fammed10[fammed10$fam%in%names(largest10),]), aes(x=distance, y=value, color=sup, group=paste(fam, variable))) + geom_hline(yintercept=0, color='gray', linetype='dashed') + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank()) + ylim(-2,2)+ xlim(0,2100)
print(ggplot(fammed10[fammed10$famsize>=10,], aes(x=distance, y=value, col=sup, group=paste(fam, variable), alpha=log10(famsize)/4)) + geom_hline(yintercept=0, color='gray', linetype='dashed') + geom_line() + scale_color_manual(values=dd.col) + facet_grid(variable~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_text(angle = 180), strip.text.x=element_blank())+ ylim(-2,2)+ xlim(0,2100))

                                                                                            
dev.off() 
