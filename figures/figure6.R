library(dplyr)
library(cowplot)
library(data.table)

source('color_palette.R')

ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-01-31.txt')
ind$sup=as.factor(ind$sup)
ind$mya=ind$age/2/3.3e-8/1e6

imp=fread('../age_model/importance_each_variable.2019-01-31.txt')
impsort=fread('../age_model/importance_combo_variable.2019-01-31.txt')

## get rmse in units of mya
imp$rmseMya
imp$rmseMya=sqrt(imp$X.IncMSE/100)/2/3.3e-8/1e6
impsort$rmseMya=sqrt(impsort$sum/100)/2/3.3e-8/1e6 ## should i be doing this before summing though?

## make a color scale for importance plots!
myColors <- brewer.pal(12,"Paired")
myColors=myColors[c(1:6,8,10,11,12)]
names(myColors) <- c('flank_base_composition', 'TE_base_composition', 'flank_methylation_mnase', 'TE_methylation_mnase',
                     'flank_closest_gene_expression', 'TE_expression', 'TE_genes', 'TE_taxonomy', 'flank_selection', 'TE_features')
colScale <- scale_colour_manual(name = "feature",values = myColors)

## and assign to categories
## assign to categories
categories=data.frame(feature=imp$feat)
categories$category=NA
categories$category[categories$feature %in% c('fam', 'sup')]='TE_taxonomy'
categories$category[categories$feature %in% c('famsize', 'pieces', 'chr', 'helLCV3p', 'helLCV5p', 'TIRlen', 'avg_ltr_length','te_bp', 'tebp', 'tespan' )]='TE_features'
categories$category[grepl('_avg_c', categories$feature) | categories$feature %in% c('shoot_prop', 'shoot_bp','root_bp', 'root_prop', 'n_shoot_hs', 'n_root_hs') ]='TE_methylation_mnase'
categories$category[categories$feature %in% c('nCG', 'nCHG', 'nCHH',  'percGC') ]='TE_base_composition'
categories$category[categories$feature %in% c('GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ENV', 'CHR', 'auton', 'pol', 'orfAA', 'helprot', 'tirprot', 'rveprot', 'tpaseprot', 'GAGfam', 'APfam', 'INTfam', 'RTfam', 'RNaseHfam', 'ENVfam', 'CHRfam', 'polfam', 'autonfam', 'helprotfam', 'rveprotfam', 'tpaseprotfam')]='TE_genes'
categories$category[grepl('_[[:digit:]]+', categories$feature) | grepl('h3k9me2_[[:digit:]]+', categories$feature) | categories$feature %in% c('flank_h3k9me2', 'flank_cg', 'flank_chg', 'flank_chh')  | (grepl('_flank', categories$feature) &  grepl('oot', categories$feature))]='flank_methylation_mnase'
categories$category[grepl('^avg_c', categories$feature)]='flank_base_composition'
categories$category[grepl('cm', categories$feature) | grepl('segsites', categories$feature) | categories$feature %in% c('closest', 'ingene', 'subgenome', 'disruptor')]='flank_selection'
categories$category[grepl('^gene_', categories$feature) ]='flank_closest_gene_expression'
categories$category[grepl('^TEfam', categories$feature) | grepl('unique', categories$feature)]='TE_expression'

## new ways to get at this - quick fix :(
categories$category[categories$feature %in% c('percGC_1kbflank','nCG_1kbflank','nCHG_1kbflank','nCHH_1kbflank','flank_bp')]='flank_base_composition'
categories$category[categories$feature %in% c('flank_n_root_hs','flank_root_bp','flank_root_prop','flank_n_shoot_hs','flank_shoot_bp','flank_shoot_prop')]='flank_methylation_mnase'
categories$category[categories$feature %in% c('segsites.bp', 'disruptor')]='TE_features'
categories[is.na(categories$category),]

                                 
library(stargazer)
stargazer(categories, summary=F, rownames=F, align=T)                                 

imp$category=mapvalues(imp$feat, from=categories$feature, to=categories$category)
#meltimp=melt(imp[rev(order(abs(imp$X.IncMSE))),][1:30,]) ## keep this managable!
#meltimpsum=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="abs") %>% summarize_if(.predicate="is.numeric", .funs="sum") %>% arrange(desc(X.IncMSE))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
#meltimpall=melt(imp[rev(order(abs(imp$X.IncMSE))),]) ## keep this managable!
meltimp=melt(imp[rev(order(abs(imp$scaled))),-'X.IncMSE'][1:30,]) ## keep this managable!
meltimpsum=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="sum", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
meltimpall=melt(imp[rev(order(abs(imp$scaled))),]) ## keep this managable!
meltimpmean=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="mean", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))


## for text:
library(stringr)
melt(imp %>% group_by(str_split_fixed(category, '_', 2)[,1]) %>% summarize_if(.predicate="is.numeric", .funs="mean", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
melt(imp %>% group_by(str_split_fixed(category, '_', 2)[,1]) %>% summarize_if(.predicate="is.numeric", .funs="sum", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))


## calc correlation with age for the 10 largest families in each sup??
feat30=data.frame(feat=as.character(imp$feat[rev(order(abs(imp$scaled)))[1:30]]))
### corr matrix by fam
for(fam in names(largest10)){
  feat30[,fam]<-NA
  for(feat in feat30$feat[-c(2,6)]){ ## have to get rid of sup and fam because they're not numeric
  feat30[feat30$feat==feat,fam]<-cor(ind$mya[ind$fam==fam], data.frame(ind)[,feat][ind$fam==fam]*1, use='na.or.complete')
    }
#  a=melt(round(cor(data.frame(ind)[ind$fam==fam, c('mya', feat30$feat)]*1, use='na.or.complete'),2))
#  
  }


## real plot
pdf(paste0('figure6.modeloutput.', Sys.Date(), '.pdf'), 30,12)

ispsc=ggplot(meltimp[meltimp$variable=='rmseMya',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimp$feat[order(meltimp[meltimp$variable=='rmseMya','value'])]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)

## make a correlation plot for alongside this one!!!! use feat30 i make above
#ispscCOR=ggplot( ###

musc=ggplot(meltimpmean[meltimpmean$variable=='rmseMya',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpmean$category[order(abs(meltimpmean[meltimpmean$variable=='rmseMya','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)

## try out some model fits!
ss=fread('../age_model/segsites.bp.2019-01-31.txt')
mSSf=ggplot(ss, aes(segsites.bp, yhat.centered/2/3.3e-8/1e6, color=sup)) + 
                geom_line(aes(group = yhat.id), alpha = 0.1, data=ss[ss$yhat.id %in% sample(1:34151, 1000),]) + 
#                stat_summary(fun.y = mean, geom = "line",  size = 0.5, aes(group=fam, color=sup), alpha=0.2)+  
                stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  
                scale_color_manual(values=dd.col)
anth=fread('../age_model/anther_avg_chh.2019-02-01.txt')
mAnthf=ggplot(anth, aes(anther_avg_chh, yhat.centered/2/3.3e-8/1e6, color=sup)) + #
                geom_line(aes(group = yhat.id), alpha = 0.1, data=anth[anth$yhat.id %in% sample(1:34151, 1000),]) + 
#                stat_summary(fun.y = mean, geom = "line",  size = 0.5, aes(group=fam, color=sup), alpha=0.2)+  
                stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  
                scale_color_manual(values=dd.col)

## plot the raw correlations
rSS=ggplot(ind, aes(x=segsites.bp, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none') + scale_color_manual(values=dd.col) + xlim(0,0.05)
rD=ggplot(ind, aes(group=disruptor, y=mya, color=sup)) + geom_boxplot() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)
rAnth=ggplot(ind, aes(x=anther_avg_chh, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)+ xlim(0,0.05)
rP=ggplot(ind, aes(x=TEfam_pollen_mature, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col) + xlim(0,2000)

r=plot_grid(rSS, mSSf + ylim(0,5), rAnth, mAnthf + ylim(0,1), labels=c('C', 'D', 'E', 'F'), ncol=2, align='v')


plot_grid(musc + theme(legend.position='none') + xlab('Reduction in square root mean squared error (Mya)') + ylab(''), 
          ispsc + theme(legend.position='none') + xlab('Reduction in square root mean squared error (Mya)') + ylab(''), 
          r, 
          ncol = 3, labels=c('A', 'B', ''), align = 'v', rel_widths=c(1,1,1.5))
          
## supp
#plot_grid(ispbysupscxM+ theme(legend.position='none'), ispbysupscxN+ theme(legend.position='none'), labels = "AUTO", ncol = 2, align = 'v', rel_widths=c(0.3,1))

dev.off()





## old stuff!
############ start plotting
##pdf(paste0('figure6.modeloutput.', Sys.Date(), '.pdf'), 30,12)
ip=ggplot(impsort, aes(y=category, x=sum)) + geom_point() + scale_y_discrete(limits=impsort$category)
#isp=ggplot(impsort[impsort$meanscaled>1.9,], aes(y=category, x=meanscaled)) + geom_point() + 
#       scale_y_discrete(limits=rev(impsort[impsort$meanscaled>1.9,] %>% arrange(desc(meanscaled)))$category)
isp=ggplot(meltimp[meltimp$variable=='X.IncMSE',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimp$feat[order(abs(meltimp[meltimp$variable=='X.IncMSE','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)
su=ggplot(meltimpsum[meltimpsum$variable=='X.IncMSE',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpsum$category[order(abs(meltimpsum[meltimpsum$variable=='X.IncMSE','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)
#su=ggplot(impsort[impsort$meanscaled>1.9,], aes(y=

ispbysup=ggplot(meltimpall[meltimpall$variable=='X.IncMSE',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimp$feat[order(meltimp[meltimp$variable=='X.IncMSE','value'])]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free')
       
       
       
       
ispsc=ggplot(meltimp[meltimp$variable=='scaled',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimp$feat[order(meltimp[meltimp$variable=='scaled','value'])]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)
susc=ggplot(meltimpsum[meltimpsum$variable=='scaled',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpsum$category[order(abs(meltimpsum[meltimpsum$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)

musc=ggplot(meltimpmean[meltimpmean$variable=='scaled',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpmean$category[order(abs(meltimpmean[meltimpmean$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)


ispbysupsc=ggplot(meltimpall[meltimpall$variable=='scaled',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free')


ispbysupscx=ggplot(meltimpall[meltimpall$variable=='scaled',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free_y')
ispbysupscxM=ggplot(meltimpall[meltimpall$variable=='scaled' & meltimpall$category=='flank_methylation_mnase',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free_y')
ispbysupscxN=ggplot(meltimpall[meltimpall$variable=='scaled' & meltimpall$category!='flank_methylation_mnase',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free_y')


plot_grid(su + theme(legend.position='none'), isp + theme(legend.position='none'), ispbysup + theme(legend.position='none'), labels = "AUTO", ncol = 3, align = 'v', rel_widths=c(1,1,2))

plot_grid(susc + theme(legend.position='none'), ispsc + theme(legend.position='none'), ispbysupsc + theme(legend.position='none'), labels = "AUTO", ncol = 3, align = 'v', rel_widths=c(1,1,2))



## plot the raw correlations
rSS=ggplot(ind, aes(x=segsites.bp, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none') + scale_color_manual(values=dd.col) + xlim(0,0.05)
rD=ggplot(ind, aes(group=disruptor, y=mya, color=sup)) + geom_boxplot() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)
rAnth=ggplot(ind, aes(x=anther_avg_chh, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)+ xlim(0,0.05)
rP=ggplot(ind, aes(x=TEfam_pollen_mature, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col) + xlim(0,2000)

r=plot_grid(rSS, rD, rAnth, rP, ncol=2, align='v')

plot_grid(musc + theme(legend.position='none'), ispsc + theme(legend.position='none'), r, ncol = 3, align = 'v', rel_widths=c(1,1,1.5))
          
## supp
plot_grid(ispbysupscxM+ theme(legend.position='none'), ispbysupscxN+ theme(legend.position='none'), labels = "AUTO", ncol = 2, align = 'v', rel_widths=c(0.3,1))

dev.off()







