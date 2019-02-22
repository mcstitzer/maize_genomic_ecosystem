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
library(ggplotify)


source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R')

## note that this expects there to be an ind data frame

ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-01-31.txt')

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


ltrgag=plot_percentagesLTR('GAG', 'GAG', invert=T)
ltrpol=plot_percentagesLTR('pol', 'Polyprotein', invert=T)
auton=plot_percentages('auton', 'All proteins', invert=T)

gagfam=plot_percentages('GAGfam', 'Family member with\nGAG', invert=T)
polfam=plot_percentages('polfam', 'Family member with\nPolyprotein', invert=T)
ltrautonfam=plot_percentagesLTR('autonfam', 'Family member with\nAll five LTR domains', invert=T)
autonfam=plot_percentages('autonfam', 'Family member with\nAll proteins', invert=T)

medianexpr=plotlargest('TEfamMedian', 'log10(Median TE \nexpression, famRPM)') + scale_y_log10()
exprpercopy=plotlargest('TEfamMedianPerCopy', 'log10(Median TE \nexpression copyRPM)') + scale_y_log10()
exprperbp=plotlargest('TEfamMedianPerBp', 'Median TE \nexpression copyRPKM') 
tetau=plotlargest('TEfam_tau', 'TE tau') + ylim(0.25,1)

geneexpr=plotlargest('gene_median', 'log10(Median expression\nclosest gene (FPKM))') + scale_y_log10()
genetau=plotlargest('gene_tau', 'Tau of closest gene') + ylim(0.25,1)
orfAA=plotlargest('orfAA', 'Longest ORF (AA)') + ylim(0,1500)
                               
## this is for the heatmap of TE expression, on a family specific, per copy basis
eer2=log2(ee[,-c(1,25:28)]+1)
eer2percopy=log2(ee[,-c(1,25:28)]/as.numeric(mapvalues(as.character(ee$fam), from=ind$fam, to=ind$famsize, warn_missing=F))+1)
rownames(eer2)=ee$fam
rownames(eer2percopy)=ee$fam
ar=data.frame(sup=substr(rownames(eer2),1,3), fam=rownames(eer2))
rownames(ar)=rownames(eer2)
hm10 = grid.grabExpr(pheatmap(eer2percopy[rownames(eer2percopy) %in% names(largest10),], treeheight_row = 0, cluster_rows=F, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[rownames(ar) %in% names(largest10),1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
hm = grid.grabExpr(pheatmap(eer2percopy[complete.cases(eer2percopy),], color=c('#000000', rev(viridis(100))), treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])

## this is for the gene heatmap, summarized at the level of TE family
ge=data.frame(ind %>% group_by(fam) %>% dplyr::summarize_at(which(grepl('gene_', colnames(ind))), funs(mean(., na.rm = TRUE))))
rownames(ge)=ge$fam
## maybe only genes within 10 kb from the te matter?
ge10kb=data.frame(ind[ind$closest<=10000,] %>% group_by(fam) %>% dplyr::summarize_at(which(grepl('gene_', colnames(ind))), funs(mean(., na.rm = TRUE))))
rownames(ge10kb)=ge10kb$fam
 
ge1kb=data.frame(ind[ind$closest<=1000,] %>% group_by(fam) %>% dplyr::summarize_at(which(grepl('gene_', colnames(ind))), funs(mean(., na.rm = TRUE))))
rownames(ge1kb)=ge1kb$fam
arg=data.frame(sup=substr(rownames(ge),1,3), fam=rownames(ge))
rownames(arg)=rownames(ge)
arg10kb=data.frame(sup=substr(rownames(ge10kb),1,3), fam=rownames(ge10kb))
rownames(arg10kb)=rownames(ge10kb)
arg1kb=data.frame(sup=substr(rownames(ge1kb),1,3), fam=rownames(ge1kb))
rownames(arg1kb)=rownames(ge1kb)

## here the number 10 stands for 10 largest fams per superfam                               
ghm10 = grid.grabExpr(pheatmap(log2(ge[rownames(ge) %in% names(largest10),2:(ncol(ge)-3)]+1), color=c('#000000', rev(viridis(100))), treeheight_row = 0, cluster_rows=F, scale='none', show_rownames=F, labels_col=gsub('gene_', '', colnames(ge)[2:(ncol(ge)-3)]), annotation_row=arg[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
ghm = grid.grabExpr(pheatmap(log2(ge[complete.cases(ge),2:(ncol(ge)-3)]+1), color=c('#000000', rev(viridis(100))), treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('gene_', '', colnames(ge)[2:(ncol(ge)-3)]), annotation_row=arg[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
#hm10 = grid.grabExpr(pheatmap(eer2[rownames(eer2) %in% names(largest10),], color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar[rownames(ar) %in% names(largest10),])[[4]])
is.grob(hm10)                               
                           
## only genes within 10 kb of the TE                               
ghm1010kb = grid.grabExpr(pheatmap(log2(ge10kb[rownames(ge10kb) %in% names(largest10),2:(ncol(ge10kb)-3)]+1), color=c('#000000', rev(viridis(100))), treeheight_row = 0, cluster_rows=F, scale='none', show_rownames=F, labels_col=gsub('gene_', '', colnames(ge10kb)[2:(ncol(ge)-3)]), annotation_row=arg10kb[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
ghm10kb = grid.grabExpr(pheatmap(log2(ge10kb[complete.cases(ge10kb),2:(ncol(ge10kb)-3)]+1), color=c('#000000', rev(viridis(100))), treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('gene_', '', colnames(ge10kb)[2:(ncol(ge)-3)]), annotation_row=arg10kb[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])

## only genes within 1 kb of the TE                               
ghm101kb = grid.grabExpr(pheatmap(log2(ge1kb[rownames(ge10kb) %in% names(largest10),2:(ncol(ge1kb)-3)]+1), color=c('#000000', rev(viridis(100))), treeheight_row = 0, cluster_rows=F, scale='none', show_rownames=F, labels_col=gsub('gene_', '', colnames(ge10kb)[2:(ncol(ge)-3)]), annotation_row=arg10kb[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
ghm1kb = grid.grabExpr(pheatmap(log2(ge1kb[complete.cases(ge1kb),2:(ncol(ge1kb)-3)]+1), color=c('#000000', rev(viridis(100))), treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('gene_', '', colnames(ge10kb)[2:(ncol(ge)-3)]), annotation_row=arg10kb[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
                             
## i can do something to scale the colors in a nicer way that puts highest turnover at right color points                       
quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(eer2percopy, n = 10)                               
hm10b = grid.grabExpr(pheatmap(eer2percopy[rownames(eer2percopy) %in% names(largest10),], breaks=mat_breaks, treeheight_row = 0, cluster_rows=F, color=c('#000000', rev(viridis(10))), scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[rownames(ar) %in% names(largest10),1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
hmb = grid.grabExpr(pheatmap(eer2percopy[complete.cases(eer2percopy),], color=c('#000000', rev(viridis(10))), breaks=mat_breaks, treeheight_row = 1, scale='none', show_rownames=F, labels_col=gsub('TEfam_', '', colnames(eer2)), annotation_row=ar[,1, drop=F], annotation_colors=list(sup=dd.col), annotation_legend=F)[[4]])
                               
                               
gagautpol=plot_grid(ltrgag, ltrauton, ltrpol, ncol=3)
#fig4 = plot_grid(auton, autonfam, gagautpol, medianexpr, tetau, labels='AUTO', ncol=1, align='v')
                               
                               
                               
                               
                               
                               
                               
                               
pdf(paste0('figure4.', Sys.Date(), '.pdf'), 16,8)
                               
fig4 = plot_grid(orfAA, auton,  gagautpol, exprpercopy, tetau, labels='AUTO', ncol=1, align='v')
                          
#hms=plot_grid(as.ggplot(hm), as.ggplot(ghm), labels='AUTO', ncol=2, align='v')
#hms10=plot_grid(as.ggplot(hm10), as.ggplot(ghm10), labels='AUTO', ncol=2, align='v')
#hms10kb=plot_grid(as.ggplot(hm), as.ggplot(ghm10kb), labels='AUTO', ncol=2, align='v')
#hms1010kb=plot_grid(as.ggplot(hm10), as.ggplot(ghm1010kb), labels='AUTO', ncol=2, align='v')
#hms1kb=plot_grid(as.ggplot(hm), as.ggplot(ghm1kb), labels='AUTO', ncol=2, align='v')
#hms101kb=plot_grid(as.ggplot(hm10), as.ggplot(ghm101kb), labels='AUTO', ncol=2, align='v')
hms3=plot_grid(as.ggplot(hm), as.ggplot(ghm), as.ggplot(ghm10kb), as.ggplot(ghm1kb), labels='AUTO', ncol=4, align='v')
hms103=plot_grid(as.ggplot(hm10), as.ggplot(ghm10), as.ggplot(ghm1010kb), as.ggplot(ghm101kb), labels='AUTO', ncol=4, align='v')

legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))


fig4hm=plot_grid(fig4, legend, as.ggplot(hm) + scale_color_manual(values=dd.col), labels=c('', '', 'G'), ncol=3, rel_widths=c(1,0.1, 0.6), align='v')
fig4hm10=plot_grid(fig4, legend, as.ggplot(hm10)+ scale_color_manual(values=dd.col), labels=c('','', 'G'), ncol=3, rel_widths=c(1,0.1, 0.6), align='v')
#fig4 ## this will have the heatmap of tissue specificity added to the right side of it
### sometimes these are weird and plot just the heatmap. I recall the pdf and replot and they're then fine. Weird!

fig4hm
fig4hm10
             
hms3
hms103

dev.off()                               


pdf(paste0('supp_genes_near_tes.', Sys.Date(), '.pdf'), 12,10)
##match genic supplement!
supp4 = plot_grid(geneexpr, genetau, labels='AUTO', ncol=1, align='v')
supp4hm=plot_grid(fig4, as.ggplot(ghm) + scale_color_manual(values=dd.col), labels=c('', 'C'), ncol=2, rel_widths=c(1, 0.6), align='v')
supp4hm10kb=plot_grid(fig4, as.ggplot(ghm10kb) + scale_color_manual(values=dd.col), labels=c('', 'C'), ncol=2, rel_widths=c(1, 0.6), align='v')

                               
                               
## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=TESUPFACTORLEVELS), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))


plots <- plot_grid(helprot, tirprot, rveprot, gag, pol, auton,  labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))
ltronly = plot_grid(ltrgag, ltrap, ltrint, ltrrt, ltrrnaseh, ltrenv, ltrchr, ltrpol, ltrauton, labels='AUTO', ncol=2, align='v')
plot_grid(ltronly, legend, ncol=2, align='v', labels='', rel_widths=c(1,0.1))
sides=plot_grid(helprot, helprotfam, tirprot, tirprotfam, rveprot, rveprotfam, gag, gagfam, pol, polfam, auton, autonfam, labels='AUTO', ncol=2, align='v')
plot_grid(sides, legend, ncol=2, align='v', labels='', rel_widths=c(1,0.1))     

                               
                               
## pdf('boxplots_of_coding.pdf')

dev.off()
                               
                               
pdf(paste0('supp_proteins.', Sys.Date(), '.pdf'), 8,14)
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

ltrchrX=plot_percentagesLTR('CHR', 'Chromodomain', invert=T, xaxis=T)

                               
plot_grid(ltrgag + ylim(0,1), ltrap + ylim(0,1), ltrint + ylim(0,1), ltrrt + ylim(0,1), ltrrnaseh + ylim(0,1), ltrenv + ylim(0,1), ltrchrX + ylim(0,1), ncol=1, align='v', labels='AUTO')
                               
dev.off()

                               
                               
                               
### hmm, autonfam is not great in the tecoding - only ltr retros
## this is a quick fix
temp=ind%>%group_by(fam) %>% summarize(autonfam=sum(auton)>0)
ind$autonfam=ind$fam %in% temp$fam[temp$autonfam]                               

ind$codingstatus=ifelse(ind$auton, 'coding copy', 'noncoding family')
ind$codingstatus[ind$codingstatus=='noncoding family' & ind$autonfam]='noncoding copy'
                               
pdf(paste0('age_v_coding.', Sys.Date(), '.pdf'), 14,8)
ggplot(ind, aes(x=autonfam, y=TEfam_tau, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=TEfam_tau, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + scale_fill_manual(values=dd.col) + facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + ylim(0,1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + ylim(0,1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_boxplot() + ylim(0,.1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_boxplot() + ylim(0,.1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)

ggplot(ind, aes(x=autonfam, y=TEfam_tau, fill=sup, group=interaction(sup, autonfam))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=TEfam_tau, fill=sup, group=interaction(sup, auton))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + scale_fill_manual(values=dd.col) + facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=autonfam, y=mya, fill=sup, group=interaction(sup, autonfam))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,.1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)
ggplot(ind, aes(x=auton, y=mya, fill=sup, group=interaction(sup, auton))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,.1)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup)

ggplot(ind, aes(x=codingstatus, y=mya, fill=sup, group=interaction(sup, codingstatus))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup) + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggplot(data.frame(ind)[ind$fam %in% names(largest10),], aes(x=auton, y=mya, fill=sup, group=interaction(fam, auton))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~fam)
ggplot(data.frame(ind)[ind$fam %in% names(largest10),], aes(x=codingstatus, y=mya, fill=sup, group=interaction(fam, codingstatus))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + ylim(0,2)+ scale_fill_manual(values=dd.col)+ facet_wrap(~fam) + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
                            
## getting weird hex errors when using coord cartesian - i'll just say in legend that the previous is artificually cut at 2 mya                        
ggplot(ind, aes(x=codingstatus, y=mya, fill=sup, group=interaction(sup, codingstatus))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + coord_cartesian(ylim=c(0,2))+ scale_fill_manual(values=dd.col)+ facet_wrap(~sup) + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
ggplot(data.frame(ind)[ind$fam %in% names(largest10),], aes(x=auton, y=mya, fill=sup, group=interaction(fam, auton))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + coord_cartesian(ylim=c(0,2))+ scale_fill_manual(values=dd.col)+ facet_wrap(~fam)
ggplot(data.frame(ind)[ind$fam %in% names(largest10),], aes(x=codingstatus, y=mya, fill=sup, group=interaction(fam, codingstatus))) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + coord_cartesian(ylim=c(0,2))+ scale_fill_manual(values=dd.col)+ facet_wrap(~fam) + theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
                             
                               
                           
                               
dev.off()

                               
                               
                               
                               
### supp expression
                               
                               
ind$gene_median_1kb=NA
ind$gene_median_1kb[ind$closest<=1000 & !is.na(ind$closest)]=ind$gene_median[ind$closest<=1000 & !is.na(ind$closest)]
                               
pdf(paste0('supp_expression.', Sys.Date(), '.pdf'), 14,12)
                             
rec=plotlargest('cmmb', 'cM/Mb')
subgenome=plot_percentages('subgenome', 'Subgenome A proportion')
geneexpr=plotlargest('gene_median', 'log10(Median expression\nclosest gene (FPKM))') + scale_y_log10()
geneexpr1kb=plotlargest('gene_median_1kb', 'log10(Median expression\nclosest gene within 1kb (FPKM))') + scale_y_log10()
genetau=plotlargest('gene_tau', 'Tau of closest gene') + ylim(0.25,1)
                               
plot_grid(rec, subgenome, geneexpr, geneexpr1kb, genetau, ncol=1, align='v', labels='AUTO')
dev.off()
                               
                               
pdf(paste0('supp_regulation.', Sys.Date(), '.pdf'), 24,12)
chg=plotlargest('nCHG', 'nCHG')
chg.flank=plotlargest('nCHG_1kbflank', 'nCHG 1kb flank')
chh=plotlargest('nCHH', 'nCHH')
chh.flank=plotlargest('nCHH_1kbflank', 'nCHH 1kb flank')
shoot=plotlargest('shoot_prop', 'Proportion shoot MNase')
shoot.flank=plotlargest('flank_shoot_prop', 'Proportion shoot MNase\nflanking')
root=plotlargest('root_prop', 'Proportion root MNase')
root.flank=plotlargest('flank_root_prop', 'Proportion root MNase\nflanking')
                               
plot_grid(chg + ylim(0,0.12), chg.flank + ylim(0,0.12), 
          chh + ylim(0,0.25), chh.flank + ylim(0,0.25), 
          shoot + ylim(0,0.4), shoot.flank + ylim(0,0.4), 
          root + ylim(0,0.4), root.flank + ylim(0,0.4), 
          ncol=2, align='v', labels='AUTO')
dev.off()                               
                               
                               
                               
                               
                               
                               
                               
