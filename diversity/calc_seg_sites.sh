#!/bin/bash -login

source ../GenomeInfo.R

mkdir -p /scratch/mstitzer
cd /scratch/mstitzer

for i in `seq 10`
do
~/software/icommands/iget /iplant/home/shared/panzea/hapmap3/hmp321/unimputed/uplifted_APGv4/hmp321_agpv4_chr${i}.vcf.gz
done


bedtools intersect -c -a $TEDISJOINED -b hmp321_agpv4_chr1.vcf.gz hmp321_agpv4_chr2.vcf.gz hmp321_agpv4_chr3.vcf.gz hmp321_agpv4_chr4.vcf.gz hmp321_agpv4_chr5.vcf.gz hmp321_agpv4_chr6.vcf.gz hmp321_agpv4_chr7.vcf.gz hmp321_agpv4_chr8.vcf.gz hmp321_agpv4_chr9.vcf.gz hmp321_agpv4_chr10.vcf.gz > ~/projects/maize_genomic_ecosystem/diversity/$GENOME.segsites.txt
bedtools intersect -c -a ~/projects/maize_genomic_ecosystem/mnase/B73.allTE.FLANK1kbEACH.gff3 -b hmp321_agpv4_chr1.vcf.gz hmp321_agpv4_chr2.vcf.gz hmp321_agpv4_chr3.vcf.gz hmp321_agpv4_chr4.vcf.gz hmp321_agpv4_chr5.vcf.gz hmp321_agpv4_chr6.vcf.gz hmp321_agpv4_chr7.vcf.gz hmp321_agpv4_chr8.vcf.gz hmp321_agpv4_chr9.vcf.gz hmp321_agpv4_chr10 > ~/projects/maize_genomic_ecosystem/diversity/$GENOME.segsites.flank.txt

rm -rf /scratch/mstitzer
