library(randomForest)
library(pdp)
library(vip)
library(cowplot)
library(data.table)
library(dplyr)
library(IPMRF)
library(reshape2)
library(plyr)
library(RColorBrewer)

source('../figures/color_palette.R')

ind=fread('B73.LTRAGE.allDescriptors.2018-12-06.txt')
ind$sup=as.factor(ind$sup)
samplerow=sample.int(n=nrow(ind), size=floor(0.5*nrow(ind)), replace=F)

## this will come in data generation soon
### Don't need anymore!!!
#ind$famsize=as.numeric(table(ind$fam)[ind$fam])
#ind=ind[,-c('segsites', 'flank_segsites')]

## make family category with <53 levels
## this works out to keeping family name for big stuff (>1053 copies, then a category of smaller fams)
#ind$famlev=factor(ind$fam, levels=c(names(tail(sort(table(ind$fam)), 52)), 'smaller'))
#ind$famlev[is.na(ind$famlev)]='smaller'

## turns out IPMRF can only deal with 32 levels (or 31?!?!?!?)
ind$famlev=factor(ind$fam, levels=c(names(tail(sort(table(ind$fam)), 30)), 'smaller'))
ind$famlev[is.na(ind$famlev)]='smaller'

## full corr matrix for interactive plot!
allcorr=melt(round(cor(ind[,c(1,5:(ncol(ind)-1))]*1, use='na.or.complete'),2))
### corr matrix by sup
suplist=vector("list", 13)
names(suplist)=names(table(ind$sup))
for(sup in names(table(ind$sup))){a=melt(round(cor(data.frame(ind)[ind$sup==sup,c(1,5:(ncol(ind)-1))]*1, use='na.or.complete'),2))
                                 suplist[[sup]]=a}
### corr matrix by fam
famlist=vector("list", 107) ## this gets you families >500
names(famlist)=head(names(rev(sort(table(ind$fam)))),107)
for(fam in head(names(rev(sort(table(ind$fam)))),107)){a=melt(round(cor(data.frame(ind)[ind$fam==fam,c(1,5:(ncol(ind)-1))]*1, use='na.or.complete'),2))
                                 famlist[[fam]]=a}
## output these so shiny can play!
for(sup in names(table(ind$sup))){allcorr[,sup]=suplist[[sup]]$value}
for(fam in head(names(rev(sort(table(ind$fam)))),107)){allcorr[,fam]=famlist[[fam]]$value}
write.table(allcorr, 'correlations_by_sup_fam.txt', row.names=F, col.names=T, sep='\t', quote=F)


## quick corr heatmap

pdf('correlation_all.pdf', 50,50)
ggplot(data=melt(round(cor(ind[,5:(ncol(ind)-1)]*1, use='na.or.complete'),2)), aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed()
## then by superfamily!
for(sup in names(table(ind$sup))){
print(ggplot(data=melt(round(cor(data.frame(ind)[ind$sup==sup,5:(ncol(ind)-1)]*1, use='na.or.complete'),2)), aes(x=Var1, y=Var2, fill=value)) + geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + ggtitle(sup))
  }
dev.off()

## split data into managable pieces!
train=ind[samplerow,]
test=ind[-samplerow,]

## remove segsites/bp
minitest=ind[sample.int(n=nrow(ind), size=34151, replace=F),]
minitest[is.na(minitest)]=-1 ## na in 
minitest.origfam=minitest$fam
minitest$fam=minitest$famlev
minitest$famlev=NULL
set.seed(1234)
subset_rf=randomForest(age~., data=minitest[,-c(2)], importance=T, keep.inbag=T, replace=F, do.trace=10, ntree=1000) ## too many factor levels for TEID
## could use strata=sup to get equal sampling by superfamily in building the tree!

minitest2=ind[sample.int(n=nrow(ind), size=500, replace=F),]
minitest2[is.na(minitest2)]=-1 ## na in 
minitest2.origfam=minitest2$fam
minitest2$fam=minitest2$famlev
minitest2$famlev=NULL


## could try drop-column importance to get at colinear groups of variables - 
##   difficult to do, needs to drop a column, retrain model, and then subtract importance value from missing model from importance from baseline model

imp=data.frame(importance(subset_rf, type=1, scale=F)) #permutation importances, see https://explained.ai/rf-importance/index.html
imp$scaled=data.frame(importance(subset_rf, type=1, scale=T))[,1]
imp$feat=rownames(imp)
imp$category=imp$feat
imp$category[grepl("flank_cg", imp$feat)]='flank_cg_methylation'
imp$category[grepl("flank_chg", imp$feat)]='flank_chg_methylation'
imp$category[grepl("flank_chh", imp$feat)]='flank_chh_methylation'
imp$category[grepl("avg_cg", imp$feat)]='te_cg_methylation'
imp$category[grepl("avg_chg", imp$feat)]='te_chg_methylation'
imp$category[grepl("avg_chh", imp$feat)]='te_chh_methylation'
impsort=imp %>% group_by(category) %>% dplyr::summarize(sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))

pdf('variable_importance_plot.pdf')
vip(subset_rf, bar=F, horizontal=F, size=1)
randomForest::varImpPlot(subset_rf)
ggplot(impsort, aes(y=category, x=sum)) + geom_point() + scale_y_discrete(limits=impsort$category)
ggplot(impsort, aes(y=category, x=meanscaled)) + geom_point() + scale_y_discrete(limits=(impsort %>% arrange(desc(meanscaled)))$category)
dev.off()


### predict on different data (like subset by superfamily or order!

## there's an is.numeric step that doesn't deal with logicals appropriately - force them to numeric here!!!
minitest3=minitest2
minitest3[,5:ncol(minitest3)]=minitest3[,5:ncol(minitest3)]*1

#newpreds=lapply(names(table(ind$sup)), function(sup) ipmrfnew(subset_rf, da=data.frame(minitest3)[minitest3$sup==sup,-c(1,2)][1:10,], ntree=1000))
#newpred=ipmrfnew(subset_rf, da=minitest3[minitest3$sup=='DTA',-c(1,2)], ntree=1000)
#sapply(1:ncol(newpred), function(x) mean(newpred[,x]))
for(sup in names(table(ind$sup))){
  newpred=ipmrfnew(subset_rf, da=data.frame(minitest3)[minitest3$sup==sup,-c(1,2)], ntree=1000) ## may need to subset to make this not over the top number of things to summarize - but some less than 10 so hard to do and I gave up with small minitest3
  temp=data.frame(meanimp=sapply(1:ncol(newpred), function(x) mean(newpred[,x])), feat=as.character(colnames(newpred)))
  imp[,sup]=as.numeric(mapvalues(imp$feat, from=temp$feat, to=temp$meanimp))
 }
#for(newpred in newpreds){
#subimp=data.frame(meanimp=sapply(1:ncol(newpred), function(x) mean(newpred[,x])), feat=as.character(colnames(newpred)))
#subimp$feat=as.character(subimp$feat)
#subimp$category=subimp$feat
#subimp$category[grepl("flank_cg", subimp$feat)]='flank_cg_methylation'
#subimp$category[grepl("flank_chg", subimp$feat)]='flank_chg_methylation'
#subimp$category[grepl("flank_chh", subimp$feat)]='flank_chh_methylation'
#subimp$category[grepl("avg_cg", subimp$feat)]='te_cg_methylation'
#subimp$category[grepl("avg_chg", subimp$feat)]='te_chg_methylation'
#subimp$category[grepl("avg_chh", subimp$feat)]='te_chh_methylation'
#subimpsort=subimp %>% group_by(category) %>% summarize(sum=sum(meanimp)) %>% arrange(desc(sum))
#print(head(subimpsort))
#}
       
## colors for poster: paired for flank vs TE
#Create a custom color scale
myColors <- brewer.pal(12,"Paired")
myColors=myColors[c(1:6,8,10,11,12)]
names(myColors) <- c('flank_base_composition', 'TE_base_composition', 'flank_methylation_mnase', 'TE_methylation_mnase',
                     'flank_closest_gene_expression', 'TE_expression', 'TE_genes', 'TE_taxonomy', 'flank_selection', 'TE_features')
colScale <- scale_colour_manual(name = "feature",values = myColors)
                                 
## assign to categories
categories=data.frame(feature=imp$feat)
categories$category=NA
categories$category[categories$feature %in% c('fam', 'sup')]='TE_taxonomy'
categories$category[categories$feature %in% c('famsize', 'pieces', 'chr', 'helLCV3p', 'helLCV5p', 'TIRlen', 'avg_ltr_length','te_bp', 'tebp', 'tespan' )]='TE_features'
categories$category[grepl('avg_', categories$feature) | categories$feature %in% c('shoot_prop', 'shoot_bp','root_bp', 'root_prop', 'n_shoot_hs', 'n_root_hs') ]='TE_methylation_mnase'
categories$category[categories$feature %in% c('nCG', 'nCHG', 'nCHH',  'percGC') ]='TE_base_composition'
categories$category[categories$feature %in% c('GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ENV', 'CHR', 'auton', 'pol', 'orfAA', 'helprot', 'tirprot', 'rveprot', 'tpaseprot')]='TE_genes'
categories$category[grepl('_[[:digit:]]+', categories$feature) | grepl('h3k9me2_[[:digit:]]+', categories$feature) | categories$feature %in% c('flank_h3k9me2', 'flank_cg', 'flank_chg', 'flank_chh')  | (grepl('_flank', categories$feature) &  grepl('oot', categories$feature))]='flank_methylation_mnase'
categories$category[grepl('avg_c', categories$feature)]='flank_base_composition'
categories$category[grepl('cm', categories$feature) | grepl('segsites', categories$feature) | categories$feature %in% c('closest', 'ingene', 'subgenome', 'disruptor')]='flank_selection'
categories$category[grepl('^gene_', categories$feature) ]='flank_closest_gene_expression'
categories$category[grepl('_combo', categories$feature) | grepl('unique', categories$feature)]='TE_expression'

## new ways to get at this - quick fix :(
categories$category[categories$feature %in% c('percGC_1kbflank','nCG_1kbflank','nCHG_1kbflank','nCHH_1kbflank','flank_bp')]='flank_base_composition'
categories$category[categories$feature %in% c('flank_n_root_hs','flank_root_bp','flank_root_prop','flank_n_shoot_hs','flank_shoot_bp','flank_shoot_prop')]='flank_methylation_mnase'
categories$category[categories$feature %in% c('segsites.bp', 'disruptor')]='TE_features'
categories[is.na(categories$category),]

imp$category=mapvalues(imp$feat, from=categories$feature, to=categories$category)
meltimp=melt(imp[rev(order(imp$X.IncMSE)),][1:20,]) ## keep this managable!
meltimpsum=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="sum") %>% arrange(desc(X.IncMSE))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))

pdf('variable_importance_bysup_plot.pdf')
ggplot(impsort, aes(y=category, x=sum)) + geom_point() + scale_y_discrete(limits=impsort$category)
ggplot(impsort, aes(y=category, x=meanscaled)) + geom_point() + scale_y_discrete(limits=(impsort %>% arrange(desc(meanscaled)))$category)

ggplot(meltimpsum[meltimpsum$variable=='X.IncMSE',], aes(x=category, y=value, fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpsum$category[order(meltimpsum[meltimpsum$variable=='X.IncMSE','value'])]) +  
       coord_flip() + scale_fill_brewer(palette='Set3') #+ colScale
ggplot(meltimpsum[meltimpsum$variable=='scaled',], aes(x=category, y=value, fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpsum$category[order(meltimpsum[meltimpsum$variable=='X.IncMSE','value'])]) +  
       coord_flip() + scale_fill_brewer(palette='Set3') #+ colScale

ggplot(meltimp[meltimp$variable!='scaled',], aes(x=feat, y=value, group=variable, fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=rev(imp$feat[order(imp$X.IncMSE)])[1:20]) + coord_flip() + scale_fill_brewer(palette='Set3') + 
       ggtitle('Feature weights for gradient boosted forest') + facet_wrap(~variable)
ggplot(meltimpsum[meltimpsum$variable!='scaled',], aes(x=category, y=value, fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=rev(wm$feature)) +  
       coord_flip() + scale_fill_brewer(palette='Set3') + 
       facet_wrap(~variable) +
       ggtitle('Feature weights for gradient boosted forest, summaries') #+                               
dev.off()                                 

### correlations based on important variables                                 
                                 
                                 
                                 
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

rm.ice=generateICE('gene_Mature_Leaf_8', subset_rf)
ggplot(rm.ice, aes(gene_Mature_Leaf_8, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(gene_Mature_Leaf_8, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(gene_Mature_Leaf_8, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(gene_Mature_Leaf_8, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)


rm.ice=generateICE('te_germk_combo', subset_rf)
ggplot(rm.ice, aes(te_germk_combo, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(te_germk_combo, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) 
ggplot(rm.ice, aes(te_germk_combo, yhat/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)
ggplot(rm.ice, aes(te_germk_combo, yhat.centered/2/3.3e-8, color=sup)) + geom_line(aes(group = yhat.id), alpha = 0.2) + stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  scale_color_manual(values=dd.col) +facet_wrap(~sup)

dev.off()













