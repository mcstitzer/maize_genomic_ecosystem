#!/bin/sh

source ../GenomeInfo.R

## note these are tabs!
awk '{print $1 "  " $4 "  " $5 "  " substr($9,4,24)}' $TEDISJOINED > $TEDISJOINED.tofasta.bed

bedtools getfasta -name -fi $GENOMEFA -bed $TEDISJOINED.tofasta.bed -fo $GENOME.disjoined.fa

## combine into one sequence per TEID
$PYTHON2 collapse_chromosomes.py $GENOME.disjoined.fa > $GENOME.eachTE.fa
