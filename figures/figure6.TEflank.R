library(RColorBrewer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(gridExtra)

source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R')
source('label_x_axis_sup.R') ## gives us grobs, and grobs.supOnly which can be used to plot an x axis with superfamily names

## note that this expects there to be an ind data frame - and can use the one generated for model building!!!
ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-10-22.txt')

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
          
largest5=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),5)))))
largest5
ind$largest5=ind$fam %in% names(largest5)                               

### get genomewide info
genomewide=fread('../base_composition/B73.basecomp_genomewide.txt')
genomewide.ss=fread('../diversity/genomewide_segsites.txt')
colnames(genomewide.ss)[1]='genomewide_bp_segsites'
genomewide.mnase=read.table(paste0('../mnase/', GENOME, '.genomewide.mnase.txt'))
colnames(genomewide.mnase)=c('n_root_hs_genomewide', 'V2', 'root_bp_genomewide', 'n_shoot_hs_genomewide', 'V5', 'shoot_bp_genomewide')


te=ind

##### prepare for decay of methylation plots
anthers=colnames(ind)[which(grepl('anther', colnames(ind)))]
SAMs=colnames(ind)[which(grepl('SAM', colnames(ind)))]
earshoots=colnames(ind)[which(grepl('earshoot', colnames(ind)))]
seedlingleafs=colnames(ind)[which(grepl('all3', colnames(ind)))]
flagleafs=colnames(ind)[which(grepl('flagleaf', colnames(ind)))]
h3k9s=colnames(ind)[which(grepl('h3k9', colnames(ind)))]

tissuecols=list(anthers, SAMs, earshoots, seedlingleafs, flagleafs, h3k9s)
names(tissuecols)=c('anther', 'SAM', 'earshoot', 'all3', 'flagleaf', 'h3k9')
tissueplots=vector("list", 6)
names(tissueplots)=names(tissuecols)
for (tissue in names(tissuecols)){
d.l=melt(data.frame(ind)[ind$fam %in% names(largest10), c(tissuecols[[tissue]], 'fam', 'sup')], id.vars=c('fam', 'sup'))
d.l$distance=as.numeric(gsub("\\D", "", d.l$variable))
if(tissue=='all3'){d.l$distance=as.numeric(gsub('^3','', d.l$distance))} ## the 3 from all3 is getting into distance!!
d.l$distance[is.na(d.l$distance)]=0
d.l$context=str_split_fixed(as.character(d.l$variable), '_', 4)[,3]
#d.l$distance[d.l$context=='h3k9me2']=as.numeric(sapply(d.l$distance[d.l$context=='h3k9me2'], function(x) substring(as.character(x), 4, nchar(as.character(x)))))
d.m=melt(data.frame(ind)[ind$famsize>=10, c(tissuecols[[tissue]], 'fam', 'sup', 'famsize')], id.vars=c('sup', 'fam', 'famsize'))
d.m$distance=as.numeric(gsub("\\D", "", d.m$variable))
if(tissue=='all3'){d.m$distance=as.numeric(gsub('^3','', d.m$distance))} ## the 3 from all3 is getting into distance!!
d.m$distance[is.na(d.m$distance)]=0
d.m$context=str_split_fixed(d.m$variable, '_', 4)[,3]
#d.m$famsize=mapvalues(d.m$fam, from=d$fam, to=d$famsize, warn_missing=F)
#d.m$famsize=as.numeric(d.m$famsize)
#d.m$distance[d.m$context=='h3k9me2']=as.numeric(sapply(d.m$distance[d.m$context=='h3k9me2'], function(x) substring(as.character(x), 4, nchar(as.character(x)))))
d.l$sup=substr(d.l$fam,1,3)
d.m$sup=substr(d.m$fam,1,3)                               
                               
d.lg=d.l %>% group_by(sup, fam, variable, distance, context) %>% dplyr::summarize(value=mean(value, na.rm=T))
d.mg=d.m%>% group_by(sup, fam, variable, distance, context, famsize) %>% dplyr::summarize(value=mean(value, na.rm=T))

tissueplots[[tissue]]=d.mg
#print(ggplot(d.mg, aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(context~sup, nrow=3) +
#                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank()) + ylab('mC proportion') + scale_alpha(guide = 'none'))
                    
}
                              
                              
                              
                              
                              
##### AND FINALLY, PLOT!                              

pdf(paste0('figure6.', Sys.Date(), '.pdf'), 32,20)
cgte.seedlingleaf=plotlargest('all3_avg_cg', 'mCG, Seedling Leaf')
chgte.seedlingleaf=plotlargest('all3_avg_chg', 'mCHG, Seedling Leaf')
chhte.seedlingleaf=plotlargest('all3_avg_chh', 'mCHH, Seeding Leaf')
cgte.anther=plotlargest('anther_avg_cg', 'mCG, Anther')
chgte.anther=plotlargest('anther_avg_chg', 'mCHG, Anther')
chhte.anther=plotlargest('anther_avg_chh', 'mCHH, Anther')
cgte.SAM=plotlargest('SAM_avg_cg', 'mCG, SAM')
chgte.SAM=plotlargest('SAM_avg_chg', 'mCHG, SAM')
chhte.SAM=plotlargest('SAM_avg_chh', 'mCHH, SAM')
cgte.earshoot=plotlargest('earshoot_avg_cg', 'mCG, Earshoot')
chgte.earshoot=plotlargest('earshoot_avg_chg', 'mCHG, Earshoot')
chhte.earshoot=plotlargest('earshoot_avg_chh', 'mCHH, Earshoot')
cgte.flagleaf=plotlargest('flagleaf_avg_cg', 'mCG, Flagleaf')
chgte.flagleaf=plotlargest('flagleaf_avg_chg', 'mCHG, Flagleaf')
chhte.flagleaf=plotlargest('flagleaf_avg_chh', 'mCHH, Flagleaf')
## these are fams with >10, colored by alpha of family size
cgflank.anther=ggplot(tissueplots[['anther']][tissueplots[['anther']]$context=='cg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCG proportion') + scale_alpha(guide = 'none') +  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chgflank.anther=ggplot(tissueplots[['anther']][tissueplots[['anther']]$context=='chg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chhflank.anther=ggplot(tissueplots[['anther']][tissueplots[['anther']]$context=='chh',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHH proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
## these are fams with >10, colored by alpha of family size
cgflank.seedlingleaf=ggplot(tissueplots[['all3']][tissueplots[['all3']]$context=='cg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chgflank.seedlingleaf=ggplot(tissueplots[['all3']][tissueplots[['all3']]$context=='chg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chhflank.seedlingleaf=ggplot(tissueplots[['all3']][tissueplots[['all3']]$context=='chh',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHH proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
## these are fams with >10, colored by alpha of family size
cgflank.SAM=ggplot(tissueplots[['SAM']][tissueplots[['SAM']]$context=='cg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chgflank.SAM=ggplot(tissueplots[['SAM']][tissueplots[['SAM']]$context=='chg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chhflank.SAM=ggplot(tissueplots[['SAM']][tissueplots[['SAM']]$context=='chh',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHH proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
## these are fams with >10, colored by alpha of family size
cgflank.earshoot=ggplot(tissueplots[['earshoot']][tissueplots[['earshoot']]$context=='cg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chgflank.earshoot=ggplot(tissueplots[['earshoot']][tissueplots[['earshoot']]$context=='chg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chhflank.earshoot=ggplot(tissueplots[['earshoot']][tissueplots[['earshoot']]$context=='chh',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHH proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
## these are fams with >10, colored by alpha of family size
cgflank.flagleaf=ggplot(tissueplots[['flagleaf']][tissueplots[['flagleaf']]$context=='cg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chgflank.flagleaf=ggplot(tissueplots[['flagleaf']][tissueplots[['flagleaf']]$context=='chg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHG proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
chhflank.flagleaf=ggplot(tissueplots[['flagleaf']][tissueplots[['flagleaf']]$context=='chh',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~factor(sup, levels=TESUPFACTORLEVELS), nrow=1) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank(), legend.position="none") + 
                              ylab('mCHH proportion') + scale_alpha(guide = 'none')+  scale_x_continuous(name='Distance from TE', breaks=c(0,1000, 2000), labels=c(0,1000,2000))
                              
gc=plotlargest('percGC', '% GC', hline=genomewide$percGC, hlinecolor='black')
cg=plotlargest('nCG', 'Proportion\nCG methylatable', hline=genomewide$nCG, hlinecolor='black')
mnase=plotlargest('shoot_prop', 'Proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen, hlinecolor='black')
diversity=plotlargest('segsites.bp', 'Proportion\nsegregating sites', hline=genomewide.ss$genomewide_bp_segsites/genomewide$seqlen, hlinecolor='black')
                              
gc.flank=plotlargest('percGC_1kbflank', 'Flanking % GC', hline=genomewide$percGC, hlinecolor='black')
cg.flank=plotlargest('nCG_1kbflank', 'Flanking proportion\nCG methylatable', hline=genomewide$nCG, hlinecolor='black')
mnase.flank=plotlargest('flank_shoot_prop', 'Flanking proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen, hlinecolor='black')
diversity.flank=plotlargest('flank_segsites.bp', 'Flanking proportion\nsegregating sites', hline=genomewide.ss$genomewide_bp_segsites/genomewide$seqlen, hlinecolor='black')

plot_grid(cgte.anther+ ylim(0,1), cgflank.anther+ ylim(0,1),
          chgte.anther+ ylim(0,1), chgflank.anther+ ylim(0,1),
          chhte.anther+ ylim(0,0.4), chhflank.anther+ ylim(0,0.4),          
          gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
#          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')

plot_grid(cgte.SAM+ ylim(0,1), cgflank.SAM+ ylim(0,1),
          chgte.SAM+ ylim(0,1), chgflank.SAM+ ylim(0,1),
          chhte.SAM+ ylim(0,0.4), chhflank.SAM+ ylim(0,0.4),          
          gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
#          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')                              
                              
plot_grid(cgte.earshoot+ ylim(0,1), cgflank.earshoot+ ylim(0,1),
          chgte.earshoot+ ylim(0,1), chgflank.earshoot+ ylim(0,1),
          chhte.earshoot+ ylim(0,0.4), chhflank.earshoot+ ylim(0,0.4),          
          gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
#          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')   
                              
 plot_grid(cgte.flagleaf+ ylim(0,1), cgflank.flagleaf+ ylim(0,1),
          chgte.flagleaf+ ylim(0,1), chgflank.flagleaf+ ylim(0,1),
          chhte.flagleaf+ ylim(0,0.4), chhflank.flagleaf+ ylim(0,0.4),          
          gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
#          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')        
                              
plot_grid(cgte.seedlingleaf+ ylim(0,1), cgflank.seedlingleaf+ ylim(0,1),
          chgte.seedlingleaf+ ylim(0,1), chgflank.seedlingleaf+ ylim(0,1),
          chhte.seedlingleaf+ ylim(0,0.4), chhflank.seedlingleaf+ ylim(0,0.4),          
          gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
#          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')                                 
dev.off()
                              
pdf(paste0('figure6.methylOnly.', Sys.Date(), '.pdf'), 34,10)
                            
plots=plot_grid(cgte.anther+ ylim(0,1), 
          chgte.anther+ ylim(0,1), 
          chhte.anther+ ylim(0,0.4),    
          ## try with these others in supplement? point was that these are 
#          gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
#          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
##          mnase+ ylim(0,0.4), mnase.flank+ ylim(0,0.4),
#          diversity+ ylim(0,0.15), diversity.flank+ ylim(0,0.15),
          labels = "AUTO", ncol = 1, align = 'v')
flankplots=plot_grid(cgflank.anther+ ylim(0,1) + ylab(''),
                     chgflank.anther+ ylim(0,1) + ylab(''), 
                     chhflank.anther+ ylim(0,0.4),      
                     ncol=1, labels=c("D", 'E', 'F'), align='v')
gg <- arrangeGrob(plots, bottom=grobs.supOnly, padding = unit(3, "line"))
ggg = arrangeGrob(flankplots, bottom=grobs.supOnly, padding=unit(3, 'line'))
#grid.newpage()
#grid.draw(gg)      
plot_grid(gg, ggg, ncol=2, align='v')
dev.off()
                              

### supplemental all tissue methylation and methylation decay!
pdf(paste0('supplemental_methylation_decay.', Sys.Date(), '.pdf'), 32,40)
                              
plot_grid(cgte.anther+ ylim(0,1), cgflank.anther+ ylim(0,1),
                    cgte.SAM+ ylim(0,1), cgflank.SAM+ ylim(0,1),
                    cgte.earshoot+ ylim(0,1), cgflank.earshoot+ ylim(0,1),
                    cgte.flagleaf+ ylim(0,1), cgflank.flagleaf+ ylim(0,1),
                    cgte.seedlingleaf+ ylim(0,1), cgflank.seedlingleaf+ ylim(0,1),
          chgte.anther+ ylim(0,1), chgflank.anther+ ylim(0,1),
                    chgte.SAM+ ylim(0,1), chgflank.SAM+ ylim(0,1),
                    chgte.earshoot+ ylim(0,1), chgflank.earshoot+ ylim(0,1),
                    chgte.flagleaf+ ylim(0,1), chgflank.flagleaf+ ylim(0,1),
                    chgte.seedlingleaf+ ylim(0,1), chgflank.seedlingleaf+ ylim(0,1),
          chhte.anther+ ylim(0,0.4), chhflank.anther+ ylim(0,0.4),   
          chhte.SAM+ ylim(0,0.4), chhflank.SAM+ ylim(0,0.4),  
          chhte.earshoot+ ylim(0,0.4), chhflank.earshoot+ ylim(0,0.4), 
          chhte.flagleaf+ ylim(0,0.4), chhflank.flagleaf+ ylim(0,0.4),
          chhte.seedlingleaf+ ylim(0,0.4), chhflank.seedlingleaf+ ylim(0,0.4),          
          labels = "AUTO", ncol = 2, align = 'v')                               
                              
dev.off()

                              
                              
## tiff format, this is S8_Fig.tif (after resizing manually in Preview to 2250 pixel width)
tiff(paste0('supplemental_methylation_decay.', Sys.Date(), '.tif'), 32,40, units='in', res=300)
                              
s8=plot_grid(cgte.anther+ ylim(0,1), cgflank.anther+ ylim(0,1),
                    cgte.SAM+ ylim(0,1), cgflank.SAM+ ylim(0,1),
                    cgte.earshoot+ ylim(0,1), cgflank.earshoot+ ylim(0,1),
                    cgte.flagleaf+ ylim(0,1), cgflank.flagleaf+ ylim(0,1),
                    cgte.seedlingleaf+ ylim(0,1), cgflank.seedlingleaf+ ylim(0,1),
          chgte.anther+ ylim(0,1), chgflank.anther+ ylim(0,1),
                    chgte.SAM+ ylim(0,1), chgflank.SAM+ ylim(0,1),
                    chgte.earshoot+ ylim(0,1), chgflank.earshoot+ ylim(0,1),
                    chgte.flagleaf+ ylim(0,1), chgflank.flagleaf+ ylim(0,1),
                    chgte.seedlingleaf+ ylim(0,1), chgflank.seedlingleaf+ ylim(0,1),
          chhte.anther+ ylim(0,0.4), chhflank.anther+ ylim(0,0.4),   
          chhte.SAM+ ylim(0,0.4), chhflank.SAM+ ylim(0,0.4),  
          chhte.earshoot+ ylim(0,0.4), chhflank.earshoot+ ylim(0,0.4), 
          chhte.flagleaf+ ylim(0,0.4), chhflank.flagleaf+ ylim(0,0.4),
          chhte.seedlingleaf+ ylim(0,0.4), chhflank.seedlingleaf+ ylim(0,0.4),          
          labels = c(LETTERS, 'AA', 'AB', 'AC', 'AD'), ncol = 2, align = 'v')                               
gg <- arrangeGrob(s8, bottom=grobs.supOnlyDoubleWidth, padding = unit(3, "line"))
grid.newpage()
grid.draw(gg)      
                              
dev.off() 
                              
## supplemental base composition
pdf(paste0('pointrange_basecomp_flank.', Sys.Date(), '.pdf'), 22,14)
gc=plotlargest('percGC', '% GC', hline=genomewide$percGC)
cg=plotlargest('nCG', 'Proportion\nCG methylatable', hline=genomewide$nCG)
chg=plotlargest('nCHG', 'Proportion\nCHG methylatable', hline=genomewide$nCHG)
chh=plotlargest('nCHH', 'Proportion\nCHH methylatable', hline=genomewide$nCHH)
mnase.r=plotlargest('root_prop', 'Proportion Mnase\nhypersensitive (root)', hline=genomewide.mnase$root_bp_genomewide/genomewide$seqlen)
mnase.s=plotlargest('shoot_prop', 'Proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen)
tg=plotlargest('nTG', 'Proportion\nTG or CA dinucleotides')
diversity=plotlargest('segsites.bp', 'Proportion\nsegregating sites', hline=genomewide.ss$genomewide_bp_segsites/genomewide$seqlen, hlinecolor='black')
                              
gc.flank=plotlargest('percGC_1kbflank', 'Flanking % GC', hline=genomewide$percGC)
cg.flank=plotlargest('nCG_1kbflank', 'Flanking proportion\nCG methylatable', hline=genomewide$nCG)
chg.flank=plotlargest('nCHG_1kbflank', 'Flanking proportion\nCHG methylatable', hline=genomewide$nCHG)
chh.flank=plotlargest('nCHH_1kbflank', 'Flanking proportion\nCHH methylatable', hline=genomewide$nCHH)
mnase.flank.r=plotlargest('flank_root_prop', 'Flanking proportion Mnase\nhypersensitive (root)', hline=genomewide.mnase$root_bp_genomewide/genomewide$seqlen)
mnase.flank.s=plotlargest('flank_shoot_prop', 'Flanking proportion Mnase\nhypersensitive (shoot)', hline=genomewide.mnase$shoot_bp_genomewide/genomewide$seqlen)
tg.flank=plotlargest('nTG_1kbflank', 'Flanking proportion\nTG or CA dinucleotides')
diversity.flank=plotlargest('flank_segsites.bp', 'Flanking proportion\nsegregating sites', hline=genomewide.ss$genomewide_bp_segsites/genomewide$seqlen, hlinecolor='black')

                              
plot_grid(gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
          chg + ylim(0,0.12), chg.flank+ ylim(0,0.12), 
          chh+ ylim(0,0.2), chh.flank+ ylim(0,0.2) , 
          tg + ylim(0.05,0.15), tg.flank + ylim(0.05,0.15),
          mnase.r+ ylim(0,0.15), mnase.flank.r + ylim(0,0.15),
          mnase.s + ylim(0,0.4), mnase.flank.s + ylim(0,0.4),
          diversity + ylim(0,0.15), diversity.flank + ylim(0,0.15),
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

                              
                              
                              
                              
## tiff format, this is S9_Fig.tif (after resizing manually in Preview to 2250 pixel width)
## uses objects from above pdf block
tiff(paste0('pointrange_basecomp_flank.', Sys.Date(), '.tif'), 22,22, units='in', res=300)
s9=plot_grid(gc + ylim(0,0.8), gc.flank+ ylim(0,0.8),  
          cg+ ylim(0,0.15),cg.flank+ ylim(0,0.15), 
          chg + ylim(0,0.12), chg.flank+ ylim(0,0.12), 
          chh+ ylim(0,0.2), chh.flank+ ylim(0,0.2) , 
          tg + ylim(0.05,0.15), tg.flank + ylim(0.05,0.15),
          mnase.r+ ylim(0,0.15), mnase.flank.r + ylim(0,0.15),
          mnase.s + ylim(0,0.4), mnase.flank.s + ylim(0,0.4),
          diversity + ylim(0,0.15), diversity.flank + ylim(0,0.15),
          labels = "AUTO", ncol = 2, align = 'v')
                              
gg <- arrangeGrob(s9, bottom=grobs.supOnlyDoubleWidth, padding = unit(3, "line"))
grid.newpage()
grid.draw(gg)      

dev.off()                              
                              
                              
### supplemental methylation
pdf(paste0('supplemental_methylation.', Sys.Date(), '.pdf'), 22,14)
for (tissue in names(tissuecols)){
cgte=plot_largest(paste0(tissue, '_avg_cg'), 'mCG')
chgte=plot_largest(paste0(tissue, '_avg_chg'), 'mCHG')
chhte=plot_largest(paste0(tissue, '_avg_chh'), 'mCHH')

title=ggdraw() + draw_label(tissue, fontface='bold')
cgflank=ggplot(tissueplots[['all3']][tissueplots[['all3']]$context=='cg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank()) + 
                              ylab('mCG proportion') + scale_alpha(guide = 'none')
chgflank=ggplot(tissueplots[['all3']][tissueplots[['all3']]$context=='chg',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank()) + 
                              ylab('mCHG proportion') + scale_alpha(guide = 'none')
chhflank=ggplot(tissueplots[['all3']][tissueplots[['all3']]$context=='chh',], aes(x=distance, y=value, col=factor(sup, levels=TESUPFACTORLEVELS), group=paste(fam, context), linetype=context, alpha=log10(famsize)/4)) + geom_line() + scale_color_manual(values=dd.col, name='superfamily') + facet_wrap(~sup) +
                               theme(strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x=element_blank()) + 
                              ylab('mCHH proportion') + scale_alpha(guide = 'none')
plot_grid(title, 
          plot_grid(cgte, cgflank,
                    chgte, chgflank,
                    chhte, chhflank,          
                    labels = "AUTO", ncol = 2, align = 'v'), 
   rel_heights=c(0.1,1))


}
dev.off()

                              
                              
                              
                              
                              
                              
## cg vs tg supp
                              

pdf('cg_vs_tg.pdf',14,8)
                              
ggplot(ind, aes(x=nCG, y=nTG)) + stat_binhex() + theme(legend.position='none') + scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F)
ggplot(ind, aes(x=nCG, y=nTG, col=sup)) + stat_binhex(col=NA) + facet_wrap(~sup) + #theme(legend.position='none') +
                              scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F)

dev.off()

                              
                              
## tiff format, this is S10_Fig.tif (after resizing manually in Preview to 2250 pixel width)
tiff(paste0('cg_vs_tg.', Sys.Date(), '.tif'), 14,8, units='in', res=300)
ggplot(ind, aes(x=nCG, y=nTG, col=sup)) + stat_binhex(col=NA) + facet_wrap(~sup) + #theme(legend.position='none') +
                              scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F) + xlab('Proportion of sites in TE in CG context') + ylab('Proportion of sites in TE in TG (or CA) context') + labs(col='superfamily')

dev.off()
                              
                              
