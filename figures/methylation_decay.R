library(data.table)
library(dplyr)
library(cowplot)
library(stringr)
library(plyr)

source('../GenomeInfo.R')
source('color_palette.R')

anther=fread('../methylation/B73.TEandFlank_methylation.anther.2018-10-24.txt')
earshoot=fread('../methylation/B73.TEandFlank_methylation.earshoot.2018-10-24.txt')
flagleaf=fread('../methylation/B73.TEandFlank_methylation.flagleaf.2018-10-25.txt')
all3
sam
h3k9me2=fread()


allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
ind=merge(anther, earshoot, all=T)
ind=merge(ind, flagleaf, all=T)
#ind=merge(ind, mnase, all=T)
ind=merge(ind, allte, all=T)
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
ind$famsize=as.numeric(table(ind$fam)[ind$fam])

anthers=colnames(ind)[which(grepl('anther', colnames(ind)))]
SAMs=colnames(ind)[which(grepl('SAM', colnames(ind)))]
earshoots=colnames(ind)[which(grepl('earshoot', colnames(ind)))]
seedlingleafs=colnames(ind)[which(grepl('all3', colnames(ind)))]
h3k9s=colnames(ind)[which(grepl('h3k9', colnames(ind)))]

d.l=melt(data.frame(ind)[ind$fam %in% names(largest10), c(anthers, 'fam')], id.vars=c('fam', 'sup'))
d.l$distance=as.numeric(gsub("\\D", "", d.l$variable))
d.l$distance[is.na(d.l$distance)]=0
d.l$context=str_split_fixed(as.character(d.l$variable), '_', 4)[,3]
#d.l$distance[d.l$context=='h3k9me2']=as.numeric(sapply(d.l$distance[d.l$context=='h3k9me2'], function(x) substring(as.character(x), 4, nchar(as.character(x)))))
d.m=melt(data.frame(ind)[ind$famsize>=10, c(anthers, 'fam', 'famsize')], id.vars=c('sup', 'fam', 'famsize'))
d.m$distance=as.numeric(gsub("\\D", "", d.m$variable))
d.m$distance[is.na(d.m$distance)]=0
d.m$context=str_split_fixed(d.m$variable, '_', 4)[,3]
#d.m$famsize=mapvalues(d.m$fam, from=d$fam, to=d$famsize, warn_missing=F)
#d.m$famsize=as.numeric(d.m$famsize)
#d.m$distance[d.m$context=='h3k9me2']=as.numeric(sapply(d.m$distance[d.m$context=='h3k9me2'], function(x) substring(as.character(x), 4, nchar(as.character(x)))))
d.l$sup=substr(d.l$fam,1,3)
d.m$sup=substr(d.m$fam,1,3)                               
                               
d.lg=d.l %>% group_by(sup, fam, variable, distance, context) %>% dplyr::summarize(value=mean(value, na.rm=T))
d.mg=d.m%>% group_by(sup, fam, variable, distance, context, famsize) %>% dplyr::summarize(value=mean(value, na.rm=T))
                               
                               
                               
pdf(paste0('methylation_decay_fig.', Sys.Date(), '.pdf'), 25,5)

#for (tissue in c('anther', 'SAM', 'earshoot', 'all3')){
tissue='anther'
## ten largest families of each sup
print(ggplot(d.lg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank())+ ggtitle(tissue))
print(ggplot(d.lg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
                               geom_point(data = df2, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10
print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank())+ ggtitle(tissue))
print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
## all families >10, alpha by famsize
print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank())+ ggtitle(tissue))
print(ggplot(d.mg, aes(x=distance, y=value, col=sup, group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col) + facet_wrap(context~sup, nrow=3, scales='free_y') +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), axis.line=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
                               geom_point(data = df3, aes(x = distance, y = value), colour = "white", alpha=0) + ggtitle(tissue))
#}
dev.off()
