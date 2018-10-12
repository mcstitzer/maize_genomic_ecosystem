#!/bin/bash

source ../GenomeInfo.R

## note these are tabs!
awk 'NF<9 {next} {print $1 "	" $4 "	" $5 "	" substr($9,4,21)}' $TEDISJOINED > $GENOME.tofasta.bed

bedtools getfasta -name -fi $GENOMEFA -bed $GENOME.tofasta.bed -fo $GENOME.disjoined.fa

## combine into one sequence per TEID
python collapse_chromosomes.py $GENOME.disjoined.fa > $GENOME.eachTE.fa
