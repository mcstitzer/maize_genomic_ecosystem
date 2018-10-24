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

ind.fam=ind %>% group_by(fam, sup) %>% summarize_if(.predicate=function(x) is.numeric(x), .funs=funs(median(., na.rm=T)))

ind.fam=merge(ind %>% group_by(fam) %>% summarize(famsize=n()), ind.fam)

pdf('famsize_vs_segsites.pdf')
ggplot(ind.fam, aes(x=famsize, y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)

ggplot(ind.fam, aes(x=log10(famsize), y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4)
ggplot(ind.fam, aes(x=log10(famsize), y=flank_segsites.bp, color=sup)) + scale_color_manual(values=dd.col) + geom_point(alpha=0.4) + facet_wrap(~sup)


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
