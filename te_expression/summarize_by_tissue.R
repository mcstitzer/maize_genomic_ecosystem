library(data.table)
library(stringr)
library(dplyr)
library(reshape2)
library(plyr)

e=fread('Walley_expresslion_18Jan19.txt')
source('../figures/color_palette.R')

colnames(e)[1]='fam'

me=melt(e, id='fam') 
me$tissue=substr(me$variable, 1, nchar(as.character(me$variable))-7)

## summarized by mean expression in each tissue
ee=dcast(me %>% group_by(fam, tissue) %>% dplyr::summarize(mean(value)), fam~tissue)
#cast(me, fam~value, mean)

colnames(ee)[2:ncol(ee)]=paste0('TEfam_', colnames(ee)[2:ncol(ee)])

## calculate TE tau, family median, per copy median
# get info about family sizes and lengths
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
techar=techar[techar$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!
techar$famsize=as.numeric(table(techar$fam)[techar$fam])
techarfam=techar %>% group_by(fam) %>% dplyr::summarize(famsize=n(), mean_tebp=mean(tebp), sum_tebp=sum(tebp))

ee$TEfamMedian=apply(ee[,2:ncol(ee)],1,median, na.rm=T) ## note this is median across means!!!!
ee$TEfamMedianPerCopy=ee$TEfamMedian/as.numeric(mapvalues(ee$fam, from=techarfam$fam, to=techarfam$famsize, warn_missing=F))
ee$TEfamMedianPerBp=ee$TEfamMedianPerCopy/as.numeric(mapvalues(ee$fam, from=techarfam$fam, to=techarfam$mean_tebp, warn_missing=F))
# and tau
tau=function(x){
    t=sum(1-x/max(x))/(length(x)-1)
  }
ee$TEfam_tau=apply(ee[,2:24], 1, tau)


write.table(ee, paste0('walley_mean_expr.', Sys.Date(), '.txt'), col.names=T, row.names=F, quote=F, sep='\t')



## add heatmap
#p9 <- as.ggplot(~image(volcano))
#p1 <- as.ggplot(~barplot(1:10)) +
#    annotate("text", x = .6, y = .5,
#             label = "Hello Base Plot", size = 5,
#             color = 'firebrick', angle=45)



pdf('heatmap_tries.pdf')
pheatmap(ee)
er=e[,2:ncol(e)]
rownames(er)=e$fam
pheatmap(er)
er2=log2(er+1)
pheatmap(er2, color=colorRampPalette(c('yellow', 'white', 'blue', 'green', 'purple'))(50))
eer2=log2(ee[,-1]+1)
rownames(eer2)=e$fam
pheatmap(eer2, color=c('#000000', c(colorRampPalette(c('yellow', 'blue', 'green', 'purple'))(50))))

pheatmap(eer2[substr(rownames(eer2),1,3)=='RLC',], color=c('#000000', c(colorRampPalette(c('yellow', 'blue', 'green', 'purple'))(50))))
pheatmap(eer2[substr(rownames(eer2),1,3)=='RLG',], color=c('#000000', c(colorRampPalette(c('yellow', 'blue', 'green', 'purple'))(50))))
pheatmap(eer2[substr(rownames(eer2),1,3)=='DHH',], color=c('#000000', c(colorRampPalette(c('yellow', 'blue', 'green', 'purple'))(50))))

pheatmap(er2, annotation_col=substr(e$fam,1,3))


#
pheatmap(eer2[substr(rownames(eer2),1,3)=='RLC',], color=c('#D3D3D3', c(colorRampPalette(c('yellow', 'blue', 'green', 'purple'))(50))))
pheatmap(eer2[substr(rownames(eer2),1,3)=='RLC',], color=c('#D3D3D3', viridis(100)))
pheatmap(eer2, color=c('#000000', rev(viridis(100))))
pheatmap(eer2[substr(rownames(eer2),1,3)=='RLC',], color=c('#D3D3D3', rev(viridis(100))))


pheatmap(eer2, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F)

pheatmap(eer2, color=c('#000000', rev(viridis(1000))), scale='none', show_rownames=F)

ar=data.frame(sup=substr(rownames(eer2),1,3))
rownames(ar)=rownames(eer2)

pheatmap(eer2, color=c('#000000', rev(viridis(10))), scale='none', show_rownames=F, annotation_row=ar)
pheatmap(eer2, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar)
pheatmap(eer2, color=c('#000000', rev(viridis(1000))), scale='none', show_rownames=F, annotation_row=ar)


pheatmap(eer2, color=c('#000000', rev(viridis(100))), scale='none', show_rownames=F, annotation_row=ar, annotation_colors=dd.col)
pheatmap(eer2, color=c('#000000', rev(viridis(1000))), scale='none', show_rownames=F, annotation_row=ar, annotation_colors=dd.col)



quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(eer2), n = 1000)
pheatmap(eer2, color=c('#000000', rev(viridis(1000))), breaks=mat_breaks, scale='none', show_rownames=F, annotation_row=ar)

dev.off()



