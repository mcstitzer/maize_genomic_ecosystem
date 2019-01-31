library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)
library(pheatmap)
library(viridis)
library(gridGraphics)
library(ggimage)

source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R')

## note that this expects there to be an ind data frame

ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-01-29.txt')

############
### THE MOST DIFFICULT PLOT EVER - set up functions
## redo largest 10 calculation
largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     

ind$mya=ind$age/3.3e-8/2/1e6                               
                               
### read in te family expression, instead of mapped to individual copies as done here:
ee=read.table('../te_expression/walley_mean_expr.2019-01-28.txt', header=T)

## all genes and stuff now in ind!
################### 


################### 
### actual plotting starts!
###################

pdf(paste0('figure4.', Sys.Date(), '.pdf'), 16,8)
helprot=plot_percentages('helprot', 'Helitron proteins', invert=T)      
tirprot=plot_percentages('tpaseprot', 'Transposase proteins', invert=T)
rveprot=plot_percentages('rveprot', 'Integrase proteins', invert=T)
anyprot=plot_percentages('anyprot', 'Protein coding', invert=T)
ltrgag=plot_percentagesLTR('GAG', 'GAG', invert=T)
gag=plot_percentages('GAG', 'GAG', invert=T)
ltrap=plot_percentagesLTR('AP', 'Aspartic Proteinase', invert=T)
ltrint=plot_percentagesLTR('INT', 'Integrase', invert=T)
ltrrt=plot_percentagesLTR('RT', 'Reverse Transcriptase', invert=T)
ltrrnaseh=plot_percentagesLTR('RNaseH', 'RNase H', invert=T)
ltrenv=plot_percentagesLTR('ENV', 'Envelope', invert=T)
ltrchr=plot_percentagesLTR('CHR', 'Chromodomain', invert=T)
ltrpol=plot_percentagesLTR('pol', 'Polyprotein', invert=T)
pol=plot_percentages('pol', 'Polyprotein', invert=T)
ltrauton=plot_percentagesLTR('auton', 'All five LTR domains', invert=T)
auton=plot_percentages('auton', 'All proteins', invert=T)
                               
                               
helprotfam=plot_percentages('helprotfam', 'Family member with\nhelitron proteins', invert=T)      
tirprotfam=plot_percentages('tpaseprotfam', 'Family member with\ntransposase proteins', invert=T)
rveprotfam=plot_percentages('rveprotfam', 'Family member with\nintegrase proteins', invert=T)
anyprotfam=plot_percentages('anyprotfam', 'Family member with\nprotein coding', invert=T)
ltrgagfam=plot_percentagesLTR('GAGfam', 'Family member with\nGAG', invert=T)
gagfam=plot_percentages('GAGfam', 'Family member with\nGAG', invert=T)
ltrapfam=plot_percentagesLTR('APfam', 'Family member with\nAspartic Proteinase', invert=T)
ltrintfam=plot_percentagesLTR('INTfam', 'Family member with\nIntegrase', invert=T)
ltrrtfam=plot_percentagesLTR('RTfam', 'Family member with\nReverse Transcriptase', invert=T)
ltrrnasehfam=plot_percentagesLTR('RNaseHfam', 'Family member with\nRNase H', invert=T)
ltrenvfam=plot_percentagesLTR('ENVfam', 'Family member with\nEnvelope', invert=T)
ltrchrfam=plot_percentagesLTR('CHRfam', 'Family member with\nChromodomain', invert=T)
ltrpolfam=plot_percentagesLTR('polfam', 'Family member with\nPolyprotein', invert=T)
polfam=plot_percentages('polfam', 'Family member with\nPolyprotein', invert=T)
ltrautonfam=plot_percentagesLTR('autonfam', 'Family member with\nAll five LTR domains', invert=T)
autonfam=plot_percentages('autonfam', 'Family member with\nAll proteins', invert=T)

medianexpr=plotlargest('TEfamMedian', 'log10(Median TE \nexpression, famRPM)') + scale_y_log10()
exprpercopy=plotlargest('TEfamMedianPerCopy', 'log10(Median TE \nexpression copyRPM)') + scale_y_log10()
exprperbp=plotlargest('TEfamMedianPerBp', 'Median TE \nexpression copyRPKM') 
tetau=plotlargest('TEfam_tau', 'TE tau') + ylim(0.25,1)

eer2=log2(ee[,-c(1,25:28)]+1)
eer2percopy=log2(ee[,-c(1,25:28)]/as.numeric(mapvalues(as.character(ee$fam), from=ind$fam, to=ind$famsize, warn_missing=F))+1)
rownames(eer2)=ee$fam
rownames(eer2percopy)=ee$fam
ar=data.frame(sup=substr(rownames(eer2),1,3), fam=rownames(eer2))
rownames(ar)=rownames(eer2)
#hm = as.ggplot(~pheatmap(eer2, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar))
#hm10 = as.ggplot(~pheatmap(eer2[rownames(eer2) %in% largest10,], color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar[rownames(ar) %in% largest10,]))
#hm.func = function() pheatmap(eer2, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar)[[4]] ## need to index the fourth!! this is where plot is stored??
#hm10.func = function() pheatmap(eer2[rownames(eer2) %in% names(largest10),], cluster_rows=F, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar[rownames(ar) %in% names(largest10),])[[4]]
#
#hm.func()
#hm10.func()
## ex. function for plotting as 
#plotfunc <- function() image(volcano) # define the function
#plotfunc() # call the function to make the plot

hm10 = grid.grabExpr(pheatmap(eer2percopy[rownames(eer2percopy) %in% names(largest10),], treeheight_row = 0, cluster_rows=F, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[rownames(ar) %in% names(largest10),1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
hm = grid.grabExpr(pheatmap(eer2percopy[complete.cases(eer2percopy),], color=c('#000000', rev(viridis(100))), treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])

#hm10 = grid.grabExpr(pheatmap(eer2[rownames(eer2) %in% names(largest10),], color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar[rownames(ar) %in% names(largest10),])[[4]])
is.grob(hm10)                               

quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(eer2percopy, n = 10)                               
hm10b = grid.grabExpr(pheatmap(eer2percopy[rownames(eer2percopy) %in% names(largest10),], breaks=mat_breaks, treeheight_row = 0, cluster_rows=F, color=c('#000000', rev(viridis(10))), scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[rownames(ar) %in% names(largest10),1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
hmb = grid.grabExpr(pheatmap(eer2percopy[complete.cases(eer2percopy),], color=c('#000000', rev(viridis(10))), breaks=mat_breaks, treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
                               
                               
gagautpol=plot_grid(ltrgag, ltrauton, ltrpol, ncol=3)
#fig4 = plot_grid(auton, autonfam, gagautpol, medianexpr, tetau, labels='AUTO', ncol=1, align='v')
fig4 = plot_grid(auton, gagautpol, exprpercopy, tetau, labels='AUTO', ncol=1, align='v')



fig4hm=plot_grid(fig4, as.ggplot(hm) + scale_color_manual(values=dd.col), labels=c('', 'E'), ncol=2, rel_widths=c(1, 0.6), align='v')
fig4hm10=plot_grid(fig4, as.ggplot(hm10)+ scale_color_manual(values=dd.col), labels=c('', 'E'), ncol=2, rel_widths=c(1, 0.6), align='v')
fig4 ## this will have the heatmap of tissue specificity added to the right side of it
fig4hm
fig4hm10
plot_grid(auton, autonfam, gagautpol, exprpercopy, tetau, labels='AUTO', ncol=1, align='v')
plot_grid(auton, autonfam, gagautpol, exprperbp, tetau, labels='AUTO', ncol=1, align='v')

## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))


plots <- plot_grid(helprot, tirprot, rveprot, gag, pol, auton,  labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))
ltronly = plot_grid(ltrgag, ltrap, ltrint, ltrrt, ltrrnaseh, ltrenv, ltrchr, ltrpol, ltrauton, labels='AUTO', ncol=2, align='v')
plot_grid(ltronly, legend, ncol=2, align='v', labels='', rel_widths=c(1,0.1))
sides=plot_grid(helprot, helprotfam, tirprot, tirprotfam, rveprot, rveprotfam, gag, gagfam, pol, polfam, auton, autonfam, labels='AUTO', ncol=2, align='v')
plot_grid(sides, legend, ncol=2, align='v', labels='', rel_widths=c(1,0.1))     

                               
                               
## pdf('boxplots_of_coding.pdf')
ggplot(ind, aes(x=autonfam, y=te_tau, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=te_tau, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + scale_fill_manual(values=dd.col) + facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + ylim(0,1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + ylim(0,1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + ylim(0,.1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + ylim(0,.1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)

                               
                               
                               

dev.off()
