library(randomForest)
library(pdp)
library(vip)
library(cowplot)
library(data.table)
library(dplyr)

source('../figures/color_palette.R')

ind=fread('B73.LTRAGE.allDescriptors.2018-11-29.txt')
ind$sup=as.factor(ind$sup)
samplerow=sample.int(n=nrow(ind), size=floor(0.5*nrow(ind)), replace=F)

## this will come in data generation soon
ind$famsize=table(ind$fam)[ind$fam]
ind=ind[,-c('segsites', 'flank_segsites')]

## make family category with <53 levels
## this works out to keeping family name for big stuff (>1053 copies, then a category of smaller fams)
ind$famlev=factor(ind$fam, levels=c(names(tail(sort(table(ind$fam)), 52)), 'smaller'))
ind$famlev[is.na(ind$famlev)]='smaller'

train=ind[samplerow,]
test=ind[-samplerow,]

## remove segsites/bp
minitest=ind[sample.int(n=nrow(ind), size=34151, replace=F),]
minitest[is.na(minitest)]=-1 ## na in 
minitest.origfam=minitest$fam
minitest$fam=minitest$famlev
minitest$famlev=NULL
set.seed(1234)
subset_rf=randomForest(age~., data=minitest[,-c(2)], importance=T, do.trace=10, ntree=1000) ## too many factor levels for TEID
## could use strata=sup to get equal sampling by superfamily in building the tree!

minitest2=ind[sample.int(n=nrow(ind), size=500, replace=F),]
minitest2[is.na(minitest2)]=-1 ## na in 
minitest2.origfam=minitest2$fam
minitest2$fam=minitest2$famlev
minitest2$famlev=NULL


## could try drop-column importance to get at colinear groups of variables - 
##   difficult to do, needs to drop a column, retrain model, and then subtract importance value from missing model from importance from baseline model

imp=data.frame(importance(subset_rf, type=1, scale=F)) #permutation importances, see https://explained.ai/rf-importance/index.html
imp$feat=rownames(imp)
imp$category=imp$feat
imp$category[grepl("flank_cg", imp$feat)]='flank_cg_methylation'
imp$category[grepl("flank_chg", imp$feat)]='flank_chg_methylation'
imp$category[grepl("flank_chh", imp$feat)]='flank_chh_methylation'
imp$category[grepl("avg_cg", imp$feat)]='te_cg_methylation'
imp$category[grepl("avg_chg", imp$feat)]='te_chg_methylation'
imp$category[grepl("avg_chh", imp$feat)]='te_chh_methylation'
impsort=imp %>% group_by(category) %>% summarize(sum=sum(X.IncMSE)) %>% arrange(desc(sum))

pdf('variable_importance_plot.pdf')
vip(subset_rf, bar=F, horizontal=F, size=1)
randomForest::varImpPlot(subset_rf)
ggplot(impsort, aes(y=category, x=sum)) + geom_point() + scale_y_discrete(limits=impsort$category)
dev.off()

########### 
## partial dependences

pdf('partial_dependences.pdf')
partialPlot(subset_rf, pred.data=minitest, x.var='segsites.bp')
partialPlot(subset_rf, pred.data=minitest, x.var='ingene')
partialPlot(subset_rf, pred.data=minitest, x.var='closest')
partialPlot(subset_rf, pred.data=minitest, x.var='sup')
partialPlot(subset_rf, pred.data=minitest, x.var='flank_segsites.bp')
partialPlot(subset_rf, pred.data=minitest, x.var='flank_root_bp')
partialPlot(subset_rf, pred.data=minitest, x.var='all3_avg_chh')
partialPlot(subset_rf, pred.data=minitest, x.var='anther_flank_cg_300')
partialPlot(subset_rf, pred.data=minitest, x.var='te_bp')
partialPlot(subset_rf, pred.data=minitest, x.var='famsize')
partialPlot(subset_rf, pred.data=minitest, x.var='shoot_bp')
partialPlot(subset_rf, pred.data=minitest, x.var='shoot_prop')
partialPlot(subset_rf, pred.data=minitest, x.var='percGC')
partialPlot(subset_rf, pred.data=minitest, x.var='disruptor')
partialPlot(subset_rf, pred.data=minitest, x.var='anther_avg_chh')
partialPlot(subset_rf, pred.data=minitest, x.var='earshoot_avg_cg')
partialPlot(subset_rf, pred.data=minitest, x.var='tespan')
partialPlot(subset_rf, pred.data=minitest, x.var='nCHH')
dev.off()

##########
## ice plots
pred.ice <- function(object, newdata) predict(object, newdata)

generateICE=function(variable, model, centeredval=0){
rm.ice <- partial(model, pred.var = variable, pred.fun = pred.ice, pred.data=minitest2)
rm.ice <- rm.ice %>%
  group_by(yhat.id) %>% # perform next operation within each yhat.id
  mutate(yhat.centered = yhat - first(yhat)) # so each curve starts at yhat = 0

rm.ice$sup=minitest$sup[rm.ice$yhat.id]
rm.ice$fam=minitest$fam
return(rm.ice)
}



pdf('ice_dependences.pdf')
rm.ice=generateICE('segsites.bp', subset_rf)
ggplot(rm.ice, aes(segsites.bp, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(segsites.bp, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(segsites.bp, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(segsites.bp, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

rm.ice=generateICE('closest', subset_rf)
ggplot(rm.ice, aes(closest, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(closest, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(closest, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(closest, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

rm.ice=generateICE('famsize', subset_rf)
ggplot(rm.ice, aes(famsize, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(famsize, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(famsize, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(famsize, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

rm.ice=generateICE('percGC', subset_rf)
ggplot(rm.ice, aes(percGC, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(percGC, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(percGC, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(percGC, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

rm.ice=generateICE('disruptor', subset_rf)
ggplot(rm.ice, aes(disruptor, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(disruptor, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(disruptor, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(disruptor, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

rm.ice=generateICE('earshoot_avg_cg', subset_rf)
ggplot(rm.ice, aes(earshoot_avg_cg, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) + xlim(0,1)
ggplot(rm.ice, aes(earshoot_avg_cg, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col)  + xlim(0,1)
ggplot(rm.ice, aes(earshoot_avg_cg, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup) + xlim(0,1)
ggplot(rm.ice, aes(earshoot_avg_cg, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup) + xlim(0,1)

rm.ice=generateICE('flank_root_bp', subset_rf)
ggplot(rm.ice, aes(flank_root_bp, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col)
ggplot(rm.ice, aes(flank_root_bp, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(flank_root_bp, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(flank_root_bp, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

rm.ice=generateICE('all3_avg_chh', subset_rf)
ggplot(rm.ice, aes(all3_avg_chh, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) + xlim(0,0.2)
ggplot(rm.ice, aes(all3_avg_chh, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) + xlim(0,0.2)
ggplot(rm.ice, aes(all3_avg_chh, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)+ xlim(0,0.2)
ggplot(rm.ice, aes(all3_avg_chh, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)+ xlim(0,0.2)



dev.off()













