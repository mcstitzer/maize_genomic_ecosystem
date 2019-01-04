library(RColorBrewer)
library(treemapify)
library(cowplot)
library(data.table)
library(plyr)
library(dplyr)

source('../GenomeInfo.R')
source('color_palette.R')
source('meanSD_functions.R')

## note that this expects there to be an ind data frame

ind=fread('../age_model/B73.LTRAGE.allDescriptors.2018-12-06.txt')

############
### THE MOST DIFFICULT PLOT EVER - set up functions
## redo largest 10 calculation
largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     


### check for TEs where a family member has the protein!

ind$auton[ind$sup=='DHH']=ind$helprot[ind$sup=='DHH'] ## make sure autonomous means what I want it to!
ind$auton[substr(ind$sup,1,2)=='DT']=ind$tpaseprot[substr(ind$sup,1,2)=='DT']
ind$auton[ind$sup %in% c('RIT', 'RIL', 'RST')]=ind$rveprot[ind$sup %in% c('RIT', 'RIL', 'RST')]

ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(GAGfam=sum(GAG)>0), by='fam')

ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(APfam=sum(AP)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(INTfam=sum(INT)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(RTfam=sum(RT)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(RNaseHfam=sum(RNaseH)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(ENVfam=sum(ENV)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(CHRfam=sum(CHR)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(polfam=sum(pol)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(autonfam=sum(auton)>0), by='fam')


ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(helprotfam=sum(helprot)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(rveprotfam=sum(rveprot)>0), by='fam')
ind=merge(ind, ind%>%group_by(fam) %>% dplyr::summarize(tpaseprotfam=sum(tpaseprot)>0), by='fam')

################### 
## get median TE expression across tissues
################### 
ind$exprMedian=apply(ind[,386:408],1,median, na.rm=T)
ind$exprMedianPerCopy=ind$exprMedian/ind$famsize
ind$exprMedianPerBp=ind$exprMedianPerCopy/ind$te_bp
# and tau
tau=function(x){
    t=sum(1-x/max(x))/(length(x)-1)
  }
ind$te_tau=apply(ind[,386:408], 1, tau)

################### 
### actual plotting starts!
###################

pdf(paste0('figure4.', Sys.Date(), '.pdf'), 16,8)
helprot=plot_percentages('helprot', 'Helitron proteins', invert=T)      
tirprot=plot_percentages('tpaseprot', 'Transposase proteins', invert=T)
rveprot=plot_percentages('rveprot', 'Integrase proteins', invert=T)
anyprot=plot_percentages('anyprot', 'Protein coding', invert=T)
ltrgag=plot_percentagesLTR('GAG', 'GAG', invert=T)
gag=plot_percentages('GAG', 'GAG', invert=T)
ltrap=plot_percentagesLTR('AP', 'Aspartic Proteinase', invert=T)
ltrint=plot_percentagesLTR('INT', 'Integrase', invert=T)
ltrrt=plot_percentagesLTR('RT', 'Reverse Transcriptase', invert=T)
ltrrnaseh=plot_percentagesLTR('RNaseH', 'RNase H', invert=T)
ltrenv=plot_percentagesLTR('ENV', 'Envelope', invert=T)
ltrchr=plot_percentagesLTR('CHR', 'Chromodomain', invert=T)
ltrpol=plot_percentagesLTR('pol', 'Polyprotein', invert=T)
pol=plot_percentages('pol', 'Polyprotein', invert=T)
ltrauton=plot_percentagesLTR('auton', 'All five LTR domains', invert=T)
auton=plot_percentages('auton', 'All proteins', invert=T)
                               
                               
helprotfam=plot_percentages('helprotfam', 'Family member with\nhelitron proteins', invert=T)      
tirprotfam=plot_percentages('tpaseprotfam', 'Family member with\ntransposase proteins', invert=T)
rveprotfam=plot_percentages('rveprotfam', 'Family member with\nintegrase proteins', invert=T)
anyprotfam=plot_percentages('anyprotfam', 'Family member with\nprotein coding', invert=T)
ltrgagfam=plot_percentagesLTR('GAGfam', 'Family member with\nGAG', invert=T)
gagfam=plot_percentages('GAGfam', 'Family member with\nGAG', invert=T)
ltrapfam=plot_percentagesLTR('APfam', 'Family member with\nAspartic Proteinase', invert=T)
ltrintfam=plot_percentagesLTR('INTfam', 'Family member with\nIntegrase', invert=T)
ltrrtfam=plot_percentagesLTR('RTfam', 'Family member with\nReverse Transcriptase', invert=T)
ltrrnasehfam=plot_percentagesLTR('RNaseHfam', 'Family member with\nRNase H', invert=T)
ltrenvfam=plot_percentagesLTR('ENVfam', 'Family member with\nEnvelope', invert=T)
ltrchrfam=plot_percentagesLTR('CHRfam', 'Family member with\nChromodomain', invert=T)
ltrpolfam=plot_percentagesLTR('polfam', 'Family member with\nPolyprotein', invert=T)
polfam=plot_percentages('polfam', 'Family member with\nPolyprotein', invert=T)
ltrautonfam=plot_percentagesLTR('autonfam', 'Family member with\nAll five LTR domains', invert=T)
autonfam=plot_percentages('autonfam', 'Family member with\nAll proteins', invert=T)

medianexpr=plotlargest('exprMedian', 'log10(Median TE \nexpression, famRPM)') + scale_y_log10()
exprpercopy=plotlargest('exprMedianPerCopy', 'log10(Median TE \nexpression copyRPM)') + scale_y_log10()
exprperbp=plotlargest('exprMedianPerBp', 'Median TE \nexpression copyRPKM') 
tetau=plotlargest('te_tau', 'TE tau')



gagautpol=plot_grid(ltrgag, ltrauton, ltrpol, ncol=3)
fig4 = plot_grid(auton, autonfam, gagautpol, medianexpr, tetau, labels='AUTO', ncol=1, align='v')

fig4 ## this will have the heatmap of tissue specificity added to the right side of it
plot_grid(auton, autonfam, gagautpol, exprpercopy, tetau, labels='AUTO', ncol=1, align='v')
plot_grid(auton, autonfam, gagautpol, exprperbp, tetau, labels='AUTO', ncol=1, align='v')

## version with a legend.
legend <- get_legend( ggplot(get_largest_quantile_backgroundbox('tebp'), aes(x=x, y=median, ymin=min, ymax=max, color=factor(sup, levels=c('DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'RLC', 'RLG', 'RLX', 'RIL', 'RIT', 'RST')), fill=sup))+ geom_pointrange(size=1)+ 
                     theme(legend.title=element_blank())+ scale_color_manual(values=dd.col))


plots <- plot_grid(helprot, tirprot, rveprot, gag, pol, auton,  labels = 'AUTO', ncol = 1, align = 'v')
plot_grid(plots,legend, ncol = 2, align = 'v', labels='', rel_widths = c(1, .1))
ltronly = plot_grid(ltrgag, ltrap, ltrint, ltrrt, ltrrnaseh, ltrenv, ltrchr, ltrpol, ltrauton, labels='AUTO', ncol=2, align='v')
plot_grid(ltronly, legend, ncol=2, align='v', labels='', rel_widths=c(1,0.1))
sides=plot_grid(helprot, helprotfam, tirprot, tirprotfam, rveprot, rveprotfam, gag, gagfam, pol, polfam, auton, autonfam, labels='AUTO', ncol=2, align='v')
plot_grid(sides, legend, ncol=2, align='v', labels='', rel_widths=c(1,0.1))     


dev.off()
