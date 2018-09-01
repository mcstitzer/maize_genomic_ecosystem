module load bedtools

## I give this 64 Gb of RAM
TEFILE=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73.structuralTEv2.2018-08-31.gff3.gz
TEDISJOINED=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73.structuralTEv2.2018-08-31.disjoined.gff3.gz

## go through each tissue that has methylation data

for TISSUE in SAM flagleaf earshoot anther all3
do
  ## combine and sort all three contexts of methylation
  paste B73_${TISSUE}.wgbs.clipOverlap_BSMAP_out.txt.100.CG.bed B73_${TISSUE}.wgbs.clipOverlap_BSMAP_out.txt.100.CHG.bed B73_${TISSUE}.wgbs.clipOverlap_BSMAP_out.txt.100.CHH.bed > ${TISSUE}_allContexts.bed
  ## remove the Chr in front of chromosome names and remove header lines
  cut -c 4- ${TISSUE}_allContexts.bed | grep -v chr > ${TISSUE}_allContexts.noChr.bed
  ## methylation of windows flanking
  bedtools closest -D a -k 40 -io -a $TEFILE -b ${TISSUE}_allContexts.noContig.noChr.bed > methylationbins_closest40.${TISSUE}_Zm00001d.bed
  ## methylation of TE itself
  bedtools intersect -wa -wb -a $TEDISJOINED -b ${TISSUE}_allContexts.noContig.noChr.bed > TE_methylation_overlap.${TISSUE}_Zm00001d.bed
done


## get h3k9 ready
cut -c 4- SRR1482372_100bpcounts_sort.chr.bed | sort -k1,1 -k2,2n > h3k9me2.bed


  ## h3k9
bedtools closest -D a -k 40 -io -a $TEFILE -b h3k9me2.bed > h3k9bins_closest40.Zm00001d.bed
bedtools intersect -wa -wb -a $TEDISJOINED -b h3k9me2.bed > TE_h3k9_overlap.Zm00001d.bed


