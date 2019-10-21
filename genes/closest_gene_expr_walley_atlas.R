library(rtracklayer)
library(stringr)
library(lme4)
library(dplyr)
library(reshape2)
library(data.table)

tg=read.table('B73_closest_gene.2019-10-21.txt', header=T, sep='\t')
## this is gene expression
#e=read.table('~/projects/b73_ecology/RawData/expression/gene_rpkm_31Aug17.txt.gz', header=T)
e=fread('walley_fpkm.txt') ## from https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v4%20RNA-SEQ%20Gene%20Atlas/walley_fpkm.txt
e[is.na(e)]=0


## unfortunately, not every TE has a nearby gene we can work with
tge=tg[tg$closestgene %in% e$geneID,]
tge[,colnames(e)]=e[e$geneID %in% tg$closestgene,][match(tg$closestgene[tg$closestgene %in% e$geneID], e$geneID[e$geneID %in% tg$closestgene]),]

## the m ending denotes using em, the mean expression across all replicates


colnames(tge)[12:34]=gsub(' ', '_', colnames(tge)[12:34])   
colnames(tge)[12:34]=gsub('/', '_', colnames(tge)[12:34]) 
colnames(tge)[12:34]=gsub('-', '_', colnames(tge)[12:34]) 
colnames(tge)[12:34]=gsub('\\(', '_', colnames(tge)[12:34]) 
colnames(tge)[12:34]=gsub(')', '_', colnames(tge)[12:34]) 
colnames(tge)[12:34]=paste0('gene_', colnames(tge)[12:34])   
                             
coef.variation <- function(x) {
    sqrt(var(x))/mean(x)
}
tge$gene_coefvar=sapply(1:nrow(tge), function(x) coef.variation(unlist(tge[x,12:34])))                             
tge$gene_median=sapply(1:nrow(tge), function(x) median(unlist(tge[x,12:34]), na.rm=T))                             
    
## calculate gene tau
tau=function(x){
    t=sum(1-x/max(x))/(length(x)-1)
  }
tge$gene_tau=apply(tge[,12:34], 1, tau)
                       
                       
tge_famsize=as.numeric(table(tge$fam)[tge$fam])
tge$famsize=tge_famsize
                       
write.table(tge, paste0('B73_closest_gene_expression.maizegdbWalley.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t')                             
                     
### now repeat for syntenic
## unfortunately, not every TE has a nearby gene we can work with
tgs=tg[tg$closestgene.syntenic %in% e$geneID,]
tgs[,colnames(e)]=e[e$geneID %in% tg$closestgene.syntenic,][match(tg$closestgene.syntenic[tg$closestgene.syntenic %in% e$geneID], e$geneID[e$geneID %in% tg$closestgene.syntenic]),]

## the m ending denotes using em, the mean expression across all replicates


colnames(tgs)[12:34]=gsub(' ', '_', colnames(tgs)[12:34])   
colnames(tgs)[12:34]=gsub('/', '_', colnames(tgs)[12:34]) 
colnames(tgs)[12:34]=gsub('-', '_', colnames(tgs)[12:34]) 
colnames(tgs)[12:34]=gsub('\\(', '_', colnames(tgs)[12:34]) 
colnames(tgs)[12:34]=gsub(')', '_', colnames(tgs)[12:34]) 
colnames(tgs)[12:34]=paste0('syntenicgene_', colnames(tgs)[12:34])   
                             
coef.variation <- function(x) {
    sqrt(var(x))/mean(x)
}
tgs$syntenicgene_coefvar=sapply(1:nrow(tgs), function(x) coef.variation(unlist(tgs[x,12:34])))                             
tgs$syntenicgene_median=sapply(1:nrow(tgs), function(x) median(unlist(tgs[x,12:34]), na.rm=T))                             
    
## calculate gene tau
tau=function(x){
    t=sum(1-x/max(x))/(length(x)-1)
  }
tgs$syntenicgene_tau=apply(tgs[,12:34], 1, tau)
                       
                       
tgs_famsize=as.numeric(table(tgs$fam)[tgs$fam])
tgs$famsize=tgs_famsize
                       
write.table(tgs, paste0('B73_closest_syntenic_gene_expression.maizegdbWalley.', Sys.Date(), '.txt'), row.names=F, col.names=T, sep='\t')                             
