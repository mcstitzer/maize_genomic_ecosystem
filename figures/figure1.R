library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)

source('../GenomeInfo.R')
source('color_palette.R')

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
gene=fread('../genes/B73_closest_gene.2018-09-20.txt')

ind=merge(techar, gene, all=T)
ind$ingene=ind$closest==0

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
                     
get_largest_quantile_backgroundbox=function(feat){
 library(plyr)
 sub2=cbind(ind[ind$largest10 & !duplicated(ind$TEID),] %>% select_(feat) , ind[ind$largest10 & !duplicated(ind$TEID),] %>% select(fam))
 sub3=sub2 %>% group_by(fam) %>% summarize_all(funs(median(., na.rm=TRUE), min=quantile(., na.rm=TRUE, 0.25), max=quantile(.,na.rm=TRUE, 0.75)))
 sub3$sup=substr(sub3$fam,1,3)
 sub3$value=NA
 sub3$plt='point'
 sub1=cbind(ind %>% select_(feat) , ind %>% select(sup))
 af=sub1 %>% group_by(sup) %>% summarize_all(funs(median(., na.rm=TRUE), min=quantile(., na.rm=TRUE, 0.25), max=quantile(.,na.rm=TRUE, 0.75)))
 af$fam=af$sup
 af$value=NA
 af$plt='suppoint'
 colnames(af)[1]='fam'
 colnames(af)[5]='sup'
 colnames(af)[2:4]=paste0(colnames(af)[2:4], '_sup')
 sub3$median_sup=as.numeric(as.character(mapvalues(sub3$sup, from=af$sup, to=af$median_sup)))
 sub3$min_sup=as.numeric(as.character(mapvalues(sub3$sup, from=af$sup, to=af$max_sup)))
 sub3$max_sup=as.numeric(as.character(mapvalues(sub3$sup, from=af$sup, to=af$min_sup)))
 sub3=sub3[match(names(largest10), sub3$fam),]
 sub3$x=1:nrow(sub3)
 sub3$px=sub3$x
 sub3$x1=sub3$x+1
 sub3$px1=sub3$x1
 sub3$x[!duplicated(sub3$sup)]=sub3$x1[!duplicated(sub3$sup)]-0.75
 sub3$x1[which(!duplicated(sub3$sup))[-1]-1]=sub3$x1[which(!duplicated(sub3$sup))[-1]-1]-0.75

# sub4=merge(sub3, af, all=T)
 return(sub3)
}

## figure out percentages
get_largest_percents_backgroundbox=function(feat, invert=FALSE){
 library(plyr)
 sub2=cbind(ind[ind$largest10 & !duplicated(ind$TEID),] %>% select_(feat) , ind[ind$largest10 & !duplicated(ind$TEID),] %>% select(fam))
 classes=names(table(sub2[,feat, with=F]))
 if(invert){classes=rev(classes)}
 print(classes)
 if (length(classes)>2){
  print("warning - make sure the first class here is one you want to compare to all others")
  }
 sub3=sub2 %>% group_by(fam) %>% dplyr::summarize_all(funs(propFirst=sum(.==classes[1], na.rm=T)/sum(!is.na(.))))
 sub3$sup=substr(sub3$fam,1,3)
 sub3$plot='family'
 sub1=cbind(ind %>% select_(feat) , ind %>% select(sup))
 af=sub1 %>% group_by(sup) %>% dplyr::summarize_all(funs(propFirst=sum(.==classes[1], na.rm=T)/sum(!is.na(.))))
 af$fam=af$sup
 af$plot='superfamily'
 sub3$supperc=as.numeric(mapvalues(sub3$sup, from=af$sup, to=af$propFirst))
 sub3=sub3[match(names(largest10), sub3$fam),]
 sub3$x=1:nrow(sub3)
 sub3$px=sub3$x
 sub3$x1=sub3$x+1
 sub3$px1=sub3$x1
 sub3$x[!duplicated(sub3$sup)]=sub3$x1[!duplicated(sub3$sup)]-0.75
 sub3$x1[which(!duplicated(sub3$sup))[-1]-1]=sub3$x1[which(!duplicated(sub3$sup))[-1]-1]-0.75

# sub4=merge(sub3, af, all=T)
 return(sub3)
}


plot_percentages=function(feat, ylab='', invert=FALSE){
 ggplot(get_largest_percents_backgroundbox(feat, invert), aes(x=px, y=propFirst, fill=sup)) + 
                     geom_point(aes(color=sup), size=2) +
                      geom_col(aes(y=supperc), alpha=0.3, width=1) + 
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) +#ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
 }


### repeat for other figure of te vs flank


plotlargest=function(feat, ylab=''){
 ggplot(get_largest_quantile_backgroundbox(feat), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup)) + 
                     geom_pointrange(fatten=4/3, size=1.5) + 
                     geom_rect(aes(xmin=x, xmax=x1, fill=sup, ymin=min_sup, ymax=max_sup), alpha=0.2, colour=NA) +
                     geom_point(aes(x=px+0.5, color=sup, y=median_sup), alpha=0.5, shape="-", size=1.5) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) + #ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
 }






################### 
### actual plotting starts!
###################

pdf(paste0('pointrange_fams_quantileWideBox_selfTE.', Sys.Date(), '.pdf'), 16,8)
tel=plotlargest('tebp', 'TE Length (bp)')
#age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
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



# figure 1                               
pdf(paste0('figure1.', Sys.Date(), '.pdf'), 20,8)
tel=plotlargest('tebp', 'TE Length (bp)')
#age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
#piece=plot_percentages('pieces', 'Proportion intact')
disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
cl=plotlargest('closest', 'Distance from \ngene (bp)')
ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)

#plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(tel, cl, ingene, piece, disr + scale_x_discrete(labels=substr(names(largest10),1,3)[!duplicated(substr(names(largest10),1,3))]),  labels = "AUTO", ncol = 1, align = 'v')
## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
#plots <- plot_grid(tel, age, cl, ingene, disr ,  labels = c('B', 'C', 'D', 'E', 'F'), ncol = 1, align = 'v')
plots <- plot_grid(tel, cl, ingene, disr ,  labels = c('B', 'C', 'D', 'E'), ncol = 1, align = 'v')
plot_grid(tempfamplot, plots,legend, ncol = 3, align = 'v', labels=c('A', '', ''), scale=c(0.96,1,1), rel_widths = c(0.35, 1, .1))                              
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
pdf('supp_TE_descriptors.pdf', 16,8)
#tel=plotlargest('seqlen', 'TE Length (bp)')
#age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
piece=plot_percentages('pieces', 'Proportion intact')
#disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
#cl=plotlargest('closest', 'Distance from \ngene (bp)')
#ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)
span=plotlargest('tespan', 'TE Span (bp)')

## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup))+ geom_pointrange(size=1)+ 
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
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('seqlen'), aes(x=x, y=median, ymin=min, ymax=max, color=sup, fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
plots <- plot_grid(anyprot, ltrgag, ltrpol, ltrauton,  labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))                              
dev.off() 
                               
