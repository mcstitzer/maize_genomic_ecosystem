library(rtracklayer)
library(stringr)
library(lme4)
library(dplyr)
library(reshape2)

tg=read.table('B73_closest_gene.2018-09-20.txt', header=T, sep='\t')
## this is gene expression
e=read.table('~/projects/b73_ecology/RawData/expression/gene_rpkm_31Aug17.txt.gz', header=T)


## use mean across reps
em=mutate(e, endo_12d = rowMeans(select(e, starts_with("endo_12d")),na.rm = TRUE), 
          peri_aleurone = rowMeans(select(e, starts_with("peri_aleurone")),na.rm = TRUE), 
          endocrown = rowMeans(select(e, starts_with("endocrown")),na.rm = TRUE), 
          symdivz = rowMeans(select(e, starts_with("symdivz")),na.rm = TRUE), 
          stodivz = rowMeans(select(e, starts_with("stodivz")),na.rm = TRUE), 
          groz = rowMeans(select(e, starts_with("groz")),na.rm = TRUE), 
          emb_20d = rowMeans(select(e, starts_with("emb_20d")),na.rm = TRUE), 
          embryo = rowMeans(select(e, starts_with("embryo")),na.rm = TRUE), 
          germk = rowMeans(select(e, starts_with("germk")),na.rm = TRUE), 
          leaf_8 = rowMeans(select(e, starts_with("leaf_8")),na.rm = TRUE), 
          earp_6m = rowMeans(select(e, starts_with("earp_6m")),na.rm = TRUE), 
          earp_2m = rowMeans(select(e, starts_with("earp_2m")),na.rm = TRUE), 
          root_m = rowMeans(select(e, starts_with("root_m")),na.rm = TRUE), 
          root_e = rowMeans(select(e, starts_with("root_e")),na.rm = TRUE), 
          root_c = rowMeans(select(e, starts_with("root_c")),na.rm = TRUE), 
          p_root = rowMeans(select(e, starts_with("p_root")),na.rm = TRUE), 
          s_root = rowMeans(select(e, starts_with("s_root")),na.rm = TRUE), 
          pollen = rowMeans(select(e, starts_with("pollen")),na.rm = TRUE), 
          silk = rowMeans(select(e, starts_with("silk")),na.rm = TRUE), 
          f_spike = rowMeans(select(e, starts_with("f_spike")),na.rm = TRUE), 
          inode_6 = rowMeans(select(e, starts_with("inode_6_")),na.rm = TRUE), 
          inode_7 = rowMeans(select(e, starts_with("inode_7_")),na.rm = TRUE), 
          meristem = rowMeans(select(e, starts_with("meristem")),na.rm = TRUE))
em=em[,69:91]

medianem=em %>% summarize_if(.predicate=function(x) is.numeric(x), .funs=funs(median(., na.rm=T)))



## unfortunately, not every TE has a nearby gene we can work with
tge=tg[tg$closestgene %in% rownames(e),]
tge[,colnames(e)]=e[rownames(e) %in% tg$closestgene,][match(tg$closestgene[tg$closestgene %in% rownames(e)], rownames(e)[rownames(e) %in% tg$closestgene]),]

## the m ending denotes using em, the mean expression across all replicates
rownames(em)=rownames(e)
tgem=tg[tg$closestgene %in% rownames(em),]
tgem[,colnames(em)]=em[rownames(em) %in% tg$closestgene,][match(tg$closestgene[tg$closestgene %in% rownames(em)], rownames(em)[rownames(em) %in% tg$closestgene]),]

tge_famsize=as.numeric(table(tge$fam)[tge$fam])
tge$famsize=tge_famsize
tgem$famsize=tge_famsize

                             
colnames(tgem)[7:29]=paste0('gene_', colnames(tgem)[7:29])   
                             
coef.variation <- function(x) {
    sqrt(var(x))/mean(x)
}
tgem$gene_coefvar=sapply(1:nrow(tgem), function(x) coef.variation(unlist(tgem[x,7:29])))                             
tgem$gene_median=sapply(1:nrow(tgem), function(x) median(unlist(tgem[x,7:29])))                             
                             
write.table(tgem, 'B73v4_closest_gene_expression.txt', row.names=F, col.names=T, sep='\t')                             
                     
