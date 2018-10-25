library(rtracklayer)
library(stringr)
library(lme4)
library(dplyr)
library(reshape2)
library(data.table)

tg=read.table('B73_closest_gene.2018-09-20.txt', header=T, sep='\t')
## this is gene expression
#e=read.table('~/projects/b73_ecology/RawData/expression/gene_rpkm_31Aug17.txt.gz', header=T)
e=fread('walley_fpkm.txt') ## from https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v4%20RNA-SEQ%20Gene%20Atlas/walley_fpkm.txt



## unfortunately, not every TE has a nearby gene we can work with
tge=tg[tg$closestgene %in% e$geneID,]
tge[,colnames(e)]=e[e$geneID %in% tg$closestgene,][match(tg$closestgene[tg$closestgene %in% e$geneID], e$geneID[e$geneID %in% tg$closestgene]),]

## the m ending denotes using em, the mean expression across all replicates

tge_famsize=as.numeric(table(tge$fam)[tge$fam])
tge$famsize=tge_famsize

colnames(tge)[8:30]=gsub(' ', '_', colnames(tge)[8:30])   
colnames(tge)[8:30]=gsub('/', '_', colnames(tge)[8:30]) 
colnames(tge)[8:30]=gsub('-', '_', colnames(tge)[8:30]) 
colnames(tge)[8:30]=gsub('\\(', '_', colnames(tge)[8:30]) 
colnames(tge)[8:30]=gsub(')', '_', colnames(tge)[8:30]) 
colnames(tge)[8:30]=paste0('gene_', colnames(tge)[8:30])   
                             
coef.variation <- function(x) {
    sqrt(var(x))/mean(x)
}
tge$gene_coefvar=sapply(1:nrow(tge), function(x) coef.variation(unlist(tge[x,8:30])))                             
tge$gene_median=sapply(1:nrow(tge), function(x) median(unlist(tge[x,8:30]), na.rm=T))                             
                             
write.table(tge, 'B73_closest_gene_expression.maizegdbWalley.txt', row.names=F, col.names=T, sep='\t')                             
                     
