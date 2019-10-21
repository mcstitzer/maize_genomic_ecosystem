library(rtracklayer)
library(stringr)
library(lme4)
library(dplyr)
library(data.table)

source('../GenomeInfo.R')

## download gene annotation if needed
subgenomegff='Zea_mays.AGPv4.40.gff3.gz'  ## hahah, this is just regular genes and I didn't change the name so it still says subgenome :(
if (!file.exists(subgenomegff)) {
    download.file('ftp://ftp.ensemblgenomes.org/pub/plants/release-40/gff3/zea_mays/Zea_mays.AGPv4.40.gff3.gz' ,'Zea_mays.AGPv4.40.gff3.gz',method="auto")
#     system('gunzip Zea_mays.AGPv4.40.gff3.gz')
}

## read in genes, filter a bit
g=import.gff3('Zea_mays.AGPv4.40.gff3.gz')
mRNA=g[g$type=='mRNA',] ## i will be able to get intronic overlaps with this, comparing to the next g
g=g[g$type %in% c('five_prime_UTR', 'exon', 'CDS', 'three_prime_UTR', 'tRNA_gene', 'lincRNA', 'miRNA'),]
g$genename=str_split_fixed(unstrsplit(g$Parent), ":", 4)[,2] ## gene name with transcript numbered
g$notransc=str_split_fixed(g$genename, "_", 2)[,1] ## gene without trancript info

## read in syntenic info
synteny='sorghum3_intell_plusteff.csv'
if (!file.exists(synteny)){
	download.file('https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/14750297/sorghum3_intell_plusteff.csv', 'sorghum3_intell_plusteff.csv', method='auto')
}
sy=fread(synteny)
sygenes=unique(c(sy$maize1_v4, sy$maize2_v4))
g$syntenic=g$notransc %in% sygenes


## read in TEs
te=import.gff3(TEDISJOINED)


## find closest gene, any orientation
closest=distanceToNearest(te, g, ignore.strand=T)

te$closest=NA
te$closest[queryHits(closest)]=mcols(closest)$distance
te$closestgene=NA
te$closestgene[queryHits(closest)]=g$notransc[subjectHits(closest)]
te$closestgenetype=NA
te$closestgenetype[queryHits(closest)]=as.character(g$type)[subjectHits(closest)]
te$closestgenesyntenic=F
te$closestgenesyntenic[queryHits(closest)]=as.logical(g$syntenic)[subjectHits(closest)]

## find closest gene, same strand as TE (if TE is unstranded, this will always retrieve the closest TE on either strand!)
closest.samestrand=distanceToNearest(te, g, ignore.strand=F)

te$closest.samestrand=NA
te$closest.samestrand[queryHits(closest.samestrand)]=mcols(closest.samestrand)$distance
te$closestgene.samestrand=NA
te$closestgene.samestrand[queryHits(closest)]=g$notransc[subjectHits(closest)]
te$closestgenetype.samestrand=NA
te$closestgenetype.samestrand[queryHits(closest)]=as.character(g$type)[subjectHits(closest)]


## find closest gene, upstream (same strand) of TE (if TE is unstranded, this will always retrieve the closest TE on either strand!)
closest.upstream=precede(te, g) ## this returns the index in g that 
closest.upstreamdist=distance(te, g[closest.upstream,])

te$closest.upstream=NA
te$closest.upstream=closest.upstreamdist
te$closestgene.upstream=NA
te$closestgene.upstream=g$notransc[closest.upstream]
te$closestgenetype.upstream=NA
te$closestgenetype.upstream=as.character(g$type)[closest.upstream]

## find closest gene, downstream (same strand) of TE (if TE is unstranded, this will always retrieve the closest TE on either strand!)
#closest.downstream=precede(te, g) ## this returns the index in g that 
#closest.downstreamdist=distance(te, g[closest.downstream,])

#te$closest.downstream=NA
#te$closest.downstream=closest.downstreamdist
#te$closestgene.downstream=NA
#te$closestgene.downstream=g$notransc[closest.downstream]
#te$closestgenetype.downstream=NA
#te$closestgenetype.downstream=as.character(g$type)[closest.downstream]



## find closest gene that is syntenic
closest.syntenic=distanceToNearest(te, g[g$syntenic,], ignore.strand=T)

te$closest.syntenic=NA
te$closest.syntenic[queryHits(closest.syntenic)]=mcols(closest.syntenic)$distance
te$closestgene.syntenic=NA
te$closestgene.syntenic[queryHits(closest.syntenic)]=g$notransc[subjectHits(closest.syntenic)]
te$closestgenetype.syntenic=NA
te$closestgenetype.syntenic[queryHits(closest.syntenic)]=as.character(g$type)[subjectHits(closest.syntenic)]


## add in fam and sup!
te$sup=substr(te$ID,1,3)
te$fam=substr(te$ID,1,8)

## summarize across disjoint ranges by finding the minimum distance
tg=as.data.frame(te[!is.na(te$closest),]) %>% group_by(sup, fam, ID) %>% dplyr::summarize(closest=min(closest), 
											  closestgene=closestgene[which.min(closest)],
											  closestgenetype=closestgenetype[which.min(closest)],
											  closestgenesyntenic=closestgenesyntenic[which.min(closest)],
											  closest.syntenic=min(closest.syntenic),
											  closestgene.syntenic=ifelse(length(which.min(closest.syntenic))==1, 
														      closestgene.syntenic[which.min(closest.syntenic)], 
														      'NA'),
											  closestgenetype.syntenic=ifelse(length(which.min(closest.syntenic))==1,
															  closestgenetype.syntenic[which.min(closest.syntenic)],
															  'NA')#,
#                                                              ### same strand
#                                                              dplyr::summarize(closest.samestrand=min(closest.samestrand),
#                                                              closestgene.samestrand=closestgene.samestrand[which.min(closest.samestrand)],
#                                                              closestgenetype.samestrand=closestgenetype.samestrand[which.min(closest.samestrand)],
#                                                              ### upstream
#                                                              dplyr::summarize(closest.upstream=min(closest.upstream),
#                                                              closestgene.upstream=closestgene.samestrand[which.min(closest.upstream)],
#                                                              closestgenetype.upstream=closestgenetype.samestrand[which.min(closest.upstream)],
#                                                              ### downstream
#                                                              dplyr::summarize(closest.downstream=min(closest.downstream),
#                                                              closestgene.downstream=closestgene.downstream[which.min(closest.downstream)],
#                                                              closestgenetype.downstream=closestgenetype.downstream[which.min(closest.downstream)]

                                                              
                                                              )
tg$closestgene.syntenic[tg$closestgene.syntenic=='NA']=NA
tg$closestgenetype.syntenic[tg$closestgenetype.syntenic=='NA']=NA

colnames(tg)[3]='TEID'
write.table(tg, paste0(GENOME, '_closest_gene.', Sys.Date(), '.txt'), quote=F, sep='\t', row.names=F, col.names=T)

## could summarize across families, not implemented here to save typing!
#tg_fam=tg %>% group_by(sup, fam) %>% summarize(mean_dist=mean(closest, na.rm=T), 
#                                               median_dist=median(closest, na.rm=T), 
#                                               n_overlaps=sum(closest==0, na.rm=T))
#write.table(tg_fam, 'B74v4_closest_gene.fam.txt', quote=F, sep='\t', row.names=F, col.names=T)

#tg_sup=tg %>% group_by(sup) %>% summarize(mean_dist=mean(closest, na.rm=T), median_dist=median(closest, na.rm=T), n_overlaps=sum(closest==0, na.rm=T))
#write.table(tg_sup, 'B74v4_closest_gene.sup.txt', quote=F, sep='\t', row.names=F, col.names=T)
