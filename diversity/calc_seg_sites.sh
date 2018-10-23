#!/bin/bash -login

source ../GenomeInfo.R

mkdir -p /scratch/mstitzer
cd /scratch/mstitzer

for i in `seq 10`
do
~/software/icommands/iget /iplant/home/shared/panzea/hapmap3/hmp321/unimputed/uplifted_APGv4/hmp321_agpv4_chr${i}.vcf.gz
done

zcat *.vcf.gz | awk '/^#/ {next} {print $1, $2-1, $2}' OFS='\t' > hmp321_agpv4_segsites.bed

bedtools intersect -c -a $TEDISJOINED -b hmp321_agpv4_segsites.bed > ~/projects/maize_genomic_ecosystem/diversity/$GENOME.segsites.txt
bedtools intersect -c -a ~/projects/maize_genomic_ecosystem/mnase/B73.allTE.FLANK1kbEACH.gff3 -b hmp321_agpv4_segsites.bed > ~/projects/maize_genomic_ecosystem/diversity/$GENOME.segsites.flank.txt

rm -rf /scratch/mstitzer
