library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R')

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
colnames(gene)[3]='TEID'
tbl=fread('../te_age/tbl_age/B73_terminalbranchlength2018-10-27.txt')
ltr=fread('../te_age/B73v4_recovered_ages.txt')
ltr$TEID=gsub('B73v4', 'Zm00001d', ltr$tename)
ltr$tename=NULL


ind=merge(techar, gene, all=T)
ind$ingene=ind$closest==0
ind=merge(ind, ltr, all.x=T, by='TEID')
ind=merge(ind, tbl, all.x=T, by=c('TEID', 'fam', 'sup'))

ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!

ind$age=ind$tbl
ind$age[!is.na(ind$k2p)]=ind$k2p[!is.na(ind$k2p)]
ind$mya=ind$age/3.3e-8/2/1e6
ind$ya=ind$age/3.3e-8/2
ind$tblmya=ind$tbl/3.3e-8/2/1e6


nrow(ind)
nrow(techar)
nrow(gene)

############
### THE MOST DIFFICULT PLOT EVER - set up functions
## redo largest 10 calculation
largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
                     





################### 
### actual plotting starts!
###################

pdf(paste0('pointrange_fams_quantileWideBox_selfTE.', Sys.Date(), '.pdf'), 16,8)
tel=plotlargest('tebp', 'TE Length (bp)')
age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
#piece=plot_percentages('pieces', 'Proportion intact')
disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
cl=plotlargest('closest', 'Distance from \ngene (bp)')
ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)

plot_grid(tel, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(tel, cl, ingene, piece, disr + scale_x_discrete(labels=substr(names(largest10),1,3)[!duplicated(substr(names(largest10),1,3))]),  labels = "AUTO", ncol = 1, align = 'v')
## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
plots <- plot_grid(tel, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plots <- plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v',  rel_widths = c(1, .1))                              
dev.off()


##### add tremapify to this figure to get real figure 1

a=ind %>% group_by(fam, sup) %>% dplyr::summarize(famcopynumber=n(), fambp=sum(tebp), avg_bp=mean(tebp))
as=a %>% group_by(sup) %>% dplyr::summarize(supcopynumber=sum(famcopynumber), supbp=sum(fambp), avg_bp_sup=weighted.mean(avg_bp, famcopynumber))
fs=merge(a, as)

famplot=ggplot(fs, aes(
  area = famcopynumber,
  fill = sup,
  subgroup = paste0(sup, '\n', round(supcopynumber, digits=0), ' copies'),
  label = paste0(sup, '\n', round(supcopynumber, digits=0), ' copies'))) +
  geom_treemap(color='gray') +
  geom_treemap_subgroup_text(
    colour = "grey10",
    place = "center",
    reflow = T) +  
   scale_fill_manual(values=dd.col) +
 geom_treemap_subgroup_border()
 
## tempplot  b/c this takes forever to plot    ?? just lacks a legend.                         
tempfamplot=ggplot(fs, aes(
  area = famcopynumber,
  fill = sup,
  subgroup = paste0(sup, '\n', round(supcopynumber, digits=0), ' copies'),
  label = paste0(sup, '\n', round(supcopynumber, digits=0), ' copies'))) +
  geom_treemap(color='gray') +
  geom_treemap_subgroup_text(
    colour = "grey10",
    place = "center",
    reflow = T) +  
   scale_fill_manual(values=dd.col) +
  theme(legend.position="none")
                               
## this one uses bp instead to contrast
famplotbp=ggplot(fs, aes(
  area = famcopynumber*avg_bp,
  fill = sup,
  subgroup = paste0(sup, '\n', round(supcopynumber*avg_bp_sup/1e6, digits=0), ' Mb'),
  label = paste0(sup, '\n', round(supcopynumber*avg_bp_sup/1e6, digits=0), ' Mb'))) +
  geom_treemap(color='gray') +
  geom_treemap_subgroup_text(
    colour = "grey10",
    place = "center",
    reflow = T) +  
   scale_fill_manual(values=dd.col)+
  theme(legend.position="none")


# figure 1                               
pdf(paste0('figure1.', Sys.Date(), '.pdf'), 24,8)
tel=plotlargest('tebp', 'TE Length (bp)')
age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
#piece=plot_percentages('pieces', 'Proportion intact')
disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
cl=plotlargest('closest', 'Distance from \ngene (bp)')
ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)

#plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(tel, cl, ingene, piece, disr + scale_x_discrete(labels=substr(names(largest10),1,3)[!duplicated(substr(names(largest10),1,3))]),  labels = "AUTO", ncol = 1, align = 'v')
## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=TESUPFACTORLEVELS)))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
#plots <- plot_grid(tel, age, cl, ingene, disr ,  labels = c('B', 'C', 'D', 'E', 'F'), ncol = 1, align = 'v')
plots <- plot_grid(tel, age, cl, ingene, disr ,  labels = c('C', 'D', 'E', 'F', 'G'), ncol = 1, align = 'v')
supplots <- plot_grid(tempfamplot, famplotbp, labels=c('A', 'B'), ncol=2, align='v', scale=0.96)
plot_grid(supplots, plots,legend, ncol = 3, align = 'v', labels=c('','', ''), scale=c(0.96,1,1), rel_widths = c(0.8, 1, .1))                              
dev.off()                             


## supp figure 1
famplotbp=ggplot(fs, aes(
  area = famcopynumber*avg_bp,
  fill = sup,
  subgroup = paste0(sup, '\n', round(supcopynumber*avg_bp_sup/1e6, digits=0), ' Mb'),
  label = paste0(sup, '\n', round(supcopynumber*avg_bp_sup/1e6, digits=0), ' Mb'))) +
  geom_treemap(color='gray') +
  geom_treemap_subgroup_text(
    colour = "grey10",
    place = "center",
    reflow = T) +  
   scale_fill_manual(values=dd.col)+
  theme(legend.position="none")
# make these giant so that I can combine in keynote
#  but also make normal sized with the bp plot
pdf('supp_famsize_copies_bp.giant.pdf', 16, 32)                       
tempfamplot        
famplotbp
dev.off()
                               
pdf('supp_famsize_copies_bp.pdf', 8, 16)                       
tempfamplot        
famplotbp
dev.off()


#### may not work below this point???
############################################                               
### supp figure 2                  
############################################
pdf(paste0('supp_TE_descriptors.', Sys.Date(), '.pdf'), 16,8)
#tel=plotlargest('seqlen', 'TE Length (bp)')
#age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
piece=plot_percentages('pieces', 'Proportion intact')
#disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
#cl=plotlargest('closest', 'Distance from \ngene (bp)')
#ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)
span=plotlargest('tespan', 'TE Span (bp)')

## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
plots <- plot_grid(piece, span, labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))                              
dev.off() 
                               
                               
                               
                               
#########################
## need to invert so we're counting  percent that have the protein! (TRUE value)                               
pdf('supp_TE_coding.pdf', 16,8)
anyprot=plot_percentages('anyprot', 'Protein coding', invert=T)
ltrgag=plot_percentages('ltrgag', 'GAG', invert=T)
ltrpol=plot_percentages('ltrpol', 'Polyprotein', invert=T)
ltrauton=plot_percentages('ltrauton', 'All five LTR domains', invert=T)
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('seqlen'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
plots <- plot_grid(anyprot, ltrgag, ltrpol, ltrauton,  labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))                              
dev.off() 
                               
