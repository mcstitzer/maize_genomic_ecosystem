library(rtracklayer)
library(stringr)
library(lme4)
library(dplyr)
library(reshape2)
library(cowplot)
library(data.table)

source('../GenomeInfo.R')
source('color_palette.R')

tgem=fread('../genes/B73v4_closest_gene_expression.txt')
ge=melt(tgem, id.vars=c("sup", 'fam', 'ID', 'closest', 'closestgene', 'closestgenetype', 'famsize'))
ge$bins=cut(ge$closest, breaks=seq(0,750000, length.out=50), include.lowest=T)


supmed=tgem %>% group_by(sup) %>% summarize_if(.predicate=function(x) is.numeric(x), .funs=funs(median(., na.rm=T)))
supge=melt(supmed, id.vars=c("sup", 'closest','famsize'))

fammed=tgem %>% group_by(sup, fam) %>% summarize_if(.predicate=function(x) is.numeric(x), .funs=funs(median(., na.rm=T)))
famge=melt(fammed, id.vars=c("sup", 'fam', 'closest','famsize'))


pdf('gene_dist_expr.pdf', 20,20)
ggplot(supge, aes(x=log10(closest), y=log10(value), color=sup)) + geom_point()  + scale_color_manual(values=dd.col)
ggplot(famge, aes(x=log10(closest), y=log10(value), color=factor(variable))) + geom_point()  + scale_color_manual(values=dd.col)
ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col)
ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable)
ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y')
ggplot(ge, aes(x=bins, y=value, color=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + facet_wrap(~variable, scales='free_y') + scale_y_log10()
                                                    
#ggplot(ge, aes(x=bins, y=value, color=sup, group=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + stat_smooth(geom='loess', se=F, alpha=0.3)
#ggplot(ge, aes(x=bins, y=value, color=sup, group=sup)) + stat_summary(fun.y="mean", geom="point")  + scale_color_manual(values=dd.col) + stat_smooth(geom='loess', se=F, alpha=0.3)

dev.off()


