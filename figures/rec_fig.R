library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R') ## note this needs an ind data frame, including a tebp column for the legend to work.

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')
colnames(gene)[3]='TEID'
rec=fread('../recombination/B73_recombination.2018-10-28.txt')
colnames(rec)[1]='TEID'

ind=merge(techar, gene, all=T)
ind$ingene=ind$closest==0
ind=merge(ind, rec, all=T)
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!


largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     



nrow(ind)
nrow(techar)
nrow(gene)

pdf(paste0('supp_TE_descriptors.', Sys.Date(), '.pdf'), 16,8)
#tel=plotlargest('seqlen', 'TE Length (bp)')
#age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
piece=plot_percentages('pieces', 'Proportion intact')
#disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
#cl=plotlargest('closest', 'Distance from \ngene (bp)')
#ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)
span=plotlargest('tespan', 'TE Span (bp)')
cm=plotlargest('cmmb', 'Recombination (cM/Mb)')


## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=TESUPFACTORLEVELS)))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
plots <- plot_grid(piece, span, cm, labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))                              
dev.off() 

