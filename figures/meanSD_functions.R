library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')
source('color_palette.R')

## note that this expects there to be an ind data frame

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


plot_percentages=function(feat, ylab='', invert=FALSE, xaxis=FALSE, angle=90){
 d = get_largest_percents_backgroundbox(feat, invert)
 if(xaxis){
 ggplot(d, aes(x=factor(px), y=propFirst, fill=sup)) + 
                     geom_point(aes(color=sup), size=2) +
                      geom_col(aes(y=supperc), alpha=0.3, width=1) + 
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) +#ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_text(hjust=1, angle=angle, size=rel(0.8))) + 
                     scale_x_discrete(labels=d$fam, breaks=factor(d$px)) +
                     ylab(ylab)
  }else{
   ggplot(d, aes(x=factor(px), y=propFirst, fill=sup)) + 
                     geom_point(aes(color=sup), size=2) +
                      geom_col(aes(y=supperc), alpha=0.3, width=1) + 
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) +#ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
                     ylab(ylab)
  }
 }


### repeat for other figure of te vs flank


plotlargest=function(feat, ylab='', xaxis=FALSE, angle=90){
 d = get_largest_quantile_backgroundbox(feat)
 if(xaxis){
 ggplot(d, aes(x=factor(x), y=median, ymin=min, ymax=max, color=sup, fill=sup)) + 
                     geom_pointrange(fatten=4/3, size=1.5) + 
                     geom_rect(aes(xmin=x, xmax=x1, fill=sup, ymin=min_sup, ymax=max_sup), alpha=0.2, colour=NA) +
                     geom_point(aes(x=px+0.5, color=sup, y=median_sup), alpha=0.5, shape="-", size=1.5) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) + #ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_text(angle=angle, hjust=1, size=rel(0.8))) + 
                     scale_x_discrete(labels=d$fam, breaks=factor(d$x)) +
                     ylab(ylab)
  }else{
   ggplot(d, aes(x=factor(x), y=median, ymin=min, ymax=max, color=sup, fill=sup)) + 
                     geom_pointrange(fatten=4/3, size=1.5) + 
                     geom_rect(aes(xmin=x, xmax=x1, fill=sup, ymin=min_sup, ymax=max_sup), alpha=0.2, colour=NA) +
                     geom_point(aes(x=px+0.5, color=sup, y=median_sup), alpha=0.5, shape="-", size=1.5) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +  
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) + #ggtitle('TE length')+ 
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
  }
 }



plot_percentagesLTR=function(feat, ylab='', invert=FALSE, xaxis=FALSE, angle=90){
d=get_largest_percents_backgroundbox(feat, invert)
if(xaxis){
 ggplot(subset(d, sup%in%c('RLC', 'RLG', 'RLX')), aes(x=factor(px), y=propFirst, fill=sup)) +
                     geom_point(aes(color=sup), size=2) +
                      geom_col(aes(y=supperc), alpha=0.3, width=1) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) +#ggtitle('TE length')+
                    theme(legend.position="none", axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_text(angle=angle, hjust=1, size=rel(0.8))) + 
                     scale_x_discrete(labels=d$fam, breaks=factor(d$x)) +
                    ylab(ylab)
}else{
 ggplot(subset(d, sup%in%c('RLC', 'RLG', 'RLX')), aes(x=factor(px), y=propFirst, fill=sup)) +
                     geom_point(aes(color=sup), size=2) + 
                      geom_col(aes(y=supperc), alpha=0.3, width=1) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) +#ggtitle('TE length')+
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
}

 }


### repeat for other figure of te vs flank


plotlargestLTR=function(feat, ylab='', xaxis=FALSE, angle=90){
d=get_largest_quantile_backgroundbox(feat)
if(xaxis){
 ggplot(subset(d, sup%in%c('RLC', 'RLG', 'RLX')), aes(x=factor(x), y=median, ymin=min, ymax=max, color=sup, fill=sup)) +
                     geom_pointrange(fatten=4/3, size=1.5) +
                     geom_rect(aes(xmin=x, xmax=x1, fill=sup, ymin=min_sup, ymax=max_sup), alpha=0.2, colour=NA) +
                     geom_point(aes(x=px+0.5, color=sup, y=median_sup), alpha=0.5, shape="-", size=1.5) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) + #ggtitle('TE length')+
                     theme(legend.position="none", axis.title.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_text(angle=angle, hjust=1, size=rel(0.8))) + 
                     scale_x_discrete(labels=d$fam, breaks=factor(d$x)) +
                     ylab(ylab)
}else{
 ggplot(subset(d, sup%in%c('RLC', 'RLG', 'RLX')), aes(x=factor(x), y=median, ymin=min, ymax=max, color=sup, fill=sup)) +
                     geom_pointrange(fatten=4/3, size=1.5) +
                     geom_rect(aes(xmin=x, xmax=x1, fill=sup, ymin=min_sup, ymax=max_sup), alpha=0.2, colour=NA) +
                     geom_point(aes(x=px+0.5, color=sup, y=median_sup), alpha=0.5, shape="-", size=1.5) +
#                     geom_ribbon(aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup), alpha = 0.3)+
#                     geom_pointrange(fatten=3, size=10, shape='-', alpha=0.4, aes(x=fam, y=median_sup, ymin=min_sup, ymax=max_sup)) +
                     scale_fill_manual(values=dd.col) +  scale_color_manual(values=dd.col) + #ggtitle('TE length')+
                     theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                     ylab(ylab)
}

 }
