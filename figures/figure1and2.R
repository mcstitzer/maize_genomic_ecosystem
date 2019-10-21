library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)
library(stargazer)
library(grid)
library(gridExtra)

source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R')
source('label_x_axis_sup.R') ## gives us grobs, and grobs.supOnly which can be used to plot an x axis with superfamily names

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-10-21.txt')
ind$mya=ind$age/3.3e-8/2/1e6
ind$ya=ind$age/3.3e-8/2
indt=fread('../age_model/B73.TBLAGE.allDescriptors.2019-10-21.txt')
ind$tblmya=indt$tbl/3.3e-8/2/1e6
indt=NULL
gc()


############
### THE MOST DIFFICULT PLOT EVER - set up functions
## redo largest 10 calculation
largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
                     

largest10df=data.frame(sup=substr(names(largest10), 1,3), fam=names(largest10), famsize=largest10)
names(largest10df)=c('Superfamily', 'Family', 'Number Copies')
stargazer(largest10df, summary=F, rownames=F, align=T)


################### 
### actual plotting starts!
###################


##### Treemapify is figure 1, pointrange is figure 2. 

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
                               
## this one uses bp instead to contrast
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

                               

supplots <- plot_grid(tempfamplot, famplotbp, labels=c('A', 'B'), ncol=2, align='v', scale=0.96)

pdf(paste0('figure1.', Sys.Date(), '.pdf'), 14,8)
supplots #sup is for superfamily, not supplement!
dev.off()

# figure 2                               
pdf(paste0('figure2.', Sys.Date(), '.pdf'), 20,8)
tel=plotlargest('tebp', 'TE Length (bp)')
age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
#piece=plot_percentages('pieces', 'Proportion intact')
disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
cl=plotlargest('closest', 'Distance from \ngene (bp)')
ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)
disrX=plot_percentages('disruptor', 'Proportion in \nanother TE', xaxis=TRUE)
autonX=plot_percentages('auton', 'Proportion copies\nautonomous', xaxis=TRUE)
ageX=plotlargest('mya', 'Age \n(million years)', xaxis=TRUE) + coord_cartesian(ylim=c(0,3))
                            
                               
#plot_grid(tel, age, cl, ingene, disr ,  labels = "AUTO", ncol = 1, align = 'v')
#plot_grid(tel, cl, ingene, piece, disr + scale_x_discrete(labels=substr(names(largest10),1,3)[!duplicated(substr(names(largest10),1,3))]),  labels = "AUTO", ncol = 1, align = 'v')
## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=TESUPFACTORLEVELS)))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
#plots <- plot_grid(tel, age, cl, ingene, disr ,  labels = c('B', 'C', 'D', 'E', 'F'), ncol = 1, align = 'v')
plots <- plot_grid(tel, cl, disr, age ,  labels=c('A', 'B', 'C', 'D'), ncol = 1, align = 'v')
#supplots <- plot_grid(tempfamplot, famplotbp, labels=c('A', 'B'), ncol=2, align='v', scale=0.96)
#plot_grid(supplots, plots,legend, ncol = 3, align = 'v', labels=c('','', ''), scale=c(0.96,1,1), rel_widths = c(0.7, 1, .1))                              
plot_grid(plots,legend, ncol = 2, align = 'v', labels=c(''), scale=c(1,1), rel_widths = c(1, .1))                              

                        

                               
#x.grob <- textGrob("Common X", 
#                   gp=gpar(fontface="bold", col="blue", fontsize=15))
#add to plot
#grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
gg <- arrangeGrob(plots, bottom=grobs, padding = unit(3, "line"))
grid.newpage()
grid.draw(gg)      
                               
dev.off()                             

##### no longer needed
## supp figure 1
#famplotbp=ggplot(fs, aes(
#  area = famcopynumber*avg_bp,
#  fill = sup,
#  subgroup = paste0(sup, '\n', round(supcopynumber*avg_bp_sup/1e6, digits=0), ' Mb'),
#  label = paste0(sup, '\n', round(supcopynumber*avg_bp_sup/1e6, digits=0), ' Mb'))) +
#  geom_treemap(color='gray') +
#  geom_treemap_subgroup_text(
#    colour = "grey10",
#    place = "center",
#    reflow = T) +  
#   scale_fill_manual(values=dd.col)+
#  theme(legend.position="none")
# make these giant so that I can combine in keynote
#  but also make normal sized with the bp plot
#pdf('supp_famsize_copies_bp.giant.pdf', 16, 32)                       
#tempfamplot        
#famplotbp
#dev.off()
#                               
#pdf('supp_famsize_copies_bp.pdf', 8, 16)                       
#tempfamplot        
#famplotbp
#dev.off()


#### may not work below this point???
############################################                               
### supp figure 2                  
############################################
pdf(paste0('supp_TE_descriptors.', Sys.Date(), '.pdf'), 16,8)
#tel=plotlargest('seqlen', 'TE Length (bp)')
#age=plotlargest('mya', 'Age \n(million years)') + coord_cartesian(ylim=c(0,3))
piece=plot_percentages('pieces', 'Proportion intact')
disr=plot_percentages('disruptor', 'Proportion in \nanother TE')
ind$disruptor.samefam=ifelse(ind$disruptor.samefam==1, F, T)
disr.samefam=plot_percentages('disruptor.samefam', 'Proportion in \na family member', invert=TRUE)
#cl=plotlargest('closest', 'Distance from \ngene (bp)')
ingene=plot_percentages('ingene', 'Proportion in \ntranscript', invert=TRUE)
span=plotlargest('tespan', 'TE Span (bp)')
syngenes=plotlargest('closest.syntenic', 'Distance from \nsyntenic gene (bp)')


## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))
plots <- plot_grid(ingene, syngenes, span, piece, disr.samefam, labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1)) 
                               
gg <- arrangeGrob(plots, bottom=grobs.supOnly, padding = unit(3, "line"))
grid.newpage()
grid.draw(gg)      
         
dev.off() 
                               
                               
                               
                               

                               
