library(dplyr)
library(cowplot)

imp=fread('../age_model/importance_each_variable.2019-01-09.txt')
impsort=fread('../age_model/importance_combo_variable.2019-01-09.txt')

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
categories$category[grepl('avg_', categories$feature) | categories$feature %in% c('shoot_prop', 'shoot_bp','root_bp', 'root_prop', 'n_shoot_hs', 'n_root_hs') ]='TE_methylation_mnase'
categories$category[categories$feature %in% c('nCG', 'nCHG', 'nCHH',  'percGC') ]='TE_base_composition'
categories$category[categories$feature %in% c('GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ENV', 'CHR', 'auton', 'pol', 'orfAA', 'helprot', 'tirprot', 'rveprot', 'tpaseprot')]='TE_genes'
categories$category[grepl('_[[:digit:]]+', categories$feature) | grepl('h3k9me2_[[:digit:]]+', categories$feature) | categories$feature %in% c('flank_h3k9me2', 'flank_cg', 'flank_chg', 'flank_chh')  | (grepl('_flank', categories$feature) &  grepl('oot', categories$feature))]='flank_methylation_mnase'
#categories$category[grepl('avg_c', categories$feature)]='flank_base_composition'
categories$category[grepl('cm', categories$feature) | grepl('segsites', categories$feature) | categories$feature %in% c('closest', 'ingene', 'subgenome', 'disruptor')]='flank_selection'
categories$category[grepl('^gene_', categories$feature) ]='flank_closest_gene_expression'
categories$category[grepl('_combo', categories$feature) | grepl('unique', categories$feature)]='TE_expression'

## new ways to get at this - quick fix :(
categories$category[categories$feature %in% c('percGC_1kbflank','nCG_1kbflank','nCHG_1kbflank','nCHH_1kbflank','flank_bp')]='flank_base_composition'
categories$category[categories$feature %in% c('flank_n_root_hs','flank_root_bp','flank_root_prop','flank_n_shoot_hs','flank_shoot_bp','flank_shoot_prop')]='flank_methylation_mnase'
categories$category[categories$feature %in% c('segsites.bp', 'disruptor')]='TE_features'
categories[is.na(categories$category),]

imp$category=mapvalues(imp$feat, from=categories$feature, to=categories$category)
meltimp=melt(imp[rev(order(abs(imp$X.IncMSE))),][1:30,]) ## keep this managable!
meltimpsum=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="abs") %>% summarize_if(.predicate="is.numeric", .funs="sum") %>% arrange(desc(X.IncMSE))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
meltimpall=melt(imp[rev(order(abs(imp$X.IncMSE))),]) ## keep this managable!


############ start plotting
pdf(paste0('figure6.modeloutput.', Sys.Date(), '.pdf'), 30,12)
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

ispbysupsc=ggplot(meltimpall[meltimpall$variable=='scaled',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free')



plot_grid(su + theme(legend.position='none'), isp + theme(legend.position='none'), ispbysup + theme(legend.position='none'), labels = "AUTO", ncol = 3, align = 'v', rel_widths=c(1,1,2))

plot_grid(susc + theme(legend.position='none'), ispsc + theme(legend.position='none'), ispbysupsc + theme(legend.position='none'), labels = "AUTO", ncol = 3, align = 'v', rel_widths=c(1,1,2))


dev.off()







