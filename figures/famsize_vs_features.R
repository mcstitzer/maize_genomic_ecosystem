library(RColorBrewer)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')
source('color_palette.R')

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
basecomp=fread('../base_composition/B73_TE_methylatable.txt')
basecomp.flank=fread('../base_composition/B73_TE_methylatable.flank.txt')
colnames(basecomp.flank)[4:7]=paste0(colnames(basecomp.flank)[4:7], '_1kbflank')
mnase=fread('../mnase/B73_TEandFlank_mnase.2018-10-23.txt')

diversity=fread('../diversity/B73.segregatingsites.TEandFlank.2018-10-23.txt')

## always bring this one in!!!
allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
ind=merge(basecomp, basecomp.flank, all=T)
ind=merge(ind, diversity, all=T)
ind=merge(ind, mnase, all=T)
ind=merge(ind, allte, all=T)
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!


## or this is easier!
ind=fread('../age_model/B73.LTRAGE.allDescriptors.2018-12-06.txt')

ind.fam=ind %>% group_by(fam, sup) %>% summarize_if(.predicate=function(x) is.numeric(x), .funs=funs(median(., na.rm=T)))

ind.fam=merge(ind %>% group_by(fam) %>% dplyr::summarize(famsize=n()), ind.fam)

                                                    
pdf('famsize_vs_segsites.pdf')
ggplot(ind.fam, aes(x=famsize, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=log10(famsize), y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=log10(famsize), y=age, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=age, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

                                                    
ggplot(ind.fam, aes(x=segsites.bp, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=segsites.bp, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=nCG, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=nCG, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)
                                                    
ggplot(ind.fam, aes(x=nCG_1kbflank, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=nCG_1kbflank, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=percGC, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=percGC, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=percGC_1kbflank, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=percGC_1kbflank, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=percGC_1kbflank, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=percGC_1kbflank, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=age, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=age, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)
ggplot(ind.fam, aes(x=age, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=age, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=age, y=flank_segsites.bp, color=sup, alpha=log10(famsize))) + scale_color_manual(values=dd.col) + geom_point()
ggplot(ind.fam, aes(x=age, y=flank_segsites.bp, color=sup, alpha=log10(famsize))) + scale_color_manual(values=dd.col) + geom_point() + facet_wrap(~sup)
ggplot(ind.fam, aes(x=age, y=segsites.bp, color=sup, alpha=log10(famsize))) + scale_color_manual(values=dd.col) + geom_point()
ggplot(ind.fam, aes(x=age, y=segsites.bp, color=sup, alpha=log10(famsize))) + scale_color_manual(values=dd.col) + geom_point() + facet_wrap(~sup)


dev.off()
                                                    
                                                    
                                                    
ind$mya=ind$age/2/3.3e-8/1e6                                          
pdf('age_vs_segsites.pdf')
ggplot(ind, aes(x=mya, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind, aes(x=mya, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)
ggplot(ind, aes(x=mya, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind, aes(x=mya, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind, aes(x=mya, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + xlim(0,2)
ggplot(ind, aes(x=mya, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)+ xlim(0,2)
ggplot(ind, aes(x=mya, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)+ xlim(0,2)
ggplot(ind, aes(x=mya, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)+ xlim(0,2)
dev.off()
pdf('age_vs_segsites.fam.pdf', 20,20)
## largest10 families
ggplot(ind[ind$largest10,], aes(x=mya, y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~fam) + stat_cor(aes(color = sup), label.x = 1, method = "pearson") + xlim(0,2)
ggplot(ind[ind$largest10,], aes(x=mya, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~fam) + stat_cor(aes(color = sup), label.x = 1, method = "pearson") + xlim(0,2)

dev.off()
                                                    
pdf('age_vs_flankmethyl.pdf', 20,20)
ggplot(ind, aes(x=mya, y=all3_avg_cg - all3_flank_cg_500, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup) + stat_cor(aes(color = sup), label.x = 1, method = "pearson") + xlim(0,2)
ggplot(ind, aes(x=flank_segsites.bp, y=all3_avg_cg - all3_flank_cg_500, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup) + stat_cor(aes(color = sup), label.x = 0.2, method = "pearson")
ggplot(ind[ind$largest10,], aes(x=mya, y=all3_avg_cg - all3_flank_cg_500, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~fam) + stat_cor(aes(color = sup), label.x = 1, method = "pearson") + xlim(0,2)
ggplot(ind[ind$largest10,], aes(x=flank_segsites.bp, y=all3_avg_cg - all3_flank_cg_500, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~fam) + stat_cor(aes(color = sup), label.x = 0.2, method = "pearson")

dev.off()

pdf('famsize_vs_everything.pdf')

ggplot(ind.fam, aes(x=log10(famsize), y=tebp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=tebp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=log10(famsize), y=percGC, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=percGC, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)
ggplot(ind.fam, aes(x=log10(famsize), y=nCG, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=nCG, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)
ggplot(ind.fam, aes(x=log10(famsize), y=nCHG, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=nCHG, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)
ggplot(ind.fam, aes(x=log10(famsize), y=nCHH, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=nCHH, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

dev.off()
