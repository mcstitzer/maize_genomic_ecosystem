## loads TEFILE, TEDISJOINED, GENOMENAME, SHORTID, GENOMEFA

## TE annotations and information about the difference between full length and disjoined gff3's can be found at https://mcstitzer.github.io/maize_TEs/

## the path to the full length TE file (including the filename)
TEFILE='/home/mstitzer/projects/maize_genomic_ecosystem/B73.2018-09-19.allTE.gff3' ## this file is available at https://github.com/mcstitzer/maize_TEs/blob/master/B73.structuralTEv2.fulllength.2018-09-19.gff3.gz

## the path to the disjoined TE file (including the filename)
TEDISJOINED='/home/mstitzer/projects/maize_genomic_ecosystem/B73.2018-09-19.allTE.disjoined.gff3' ## this file is available at https://github.com/mcstitzer/maize_TEs/blob/master/B73.structuralTEv2.disjoined.2018-09-19.gff3.gz

## names of the genome and its short ID assigned by the maizegdb nomenclature committee (https://www.maizegdb.org/nomenclature)
GENOME='B73'
SHORTID='Zm00001d'

## the path to the genome fasta (including the filename)
GENOMEFA='/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73V4.both_pseudo_AND_unplaced.fa'

#### Note that this project was initated during the assembly of the B73v4 genome, so this exact fasta is not publicly available. The underlying sequence is identical to the released version, B73v4 aka B73 AGPv4.
# My fasta file closely matches the maizegdb fasta available here:
# https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0/Zm-B73-REFERENCE-GRAMENE-4.0.fa.gz
# A minor difference is that chromosomes are named Chr1-Chr10, while in the file I used, chromosomes were named 1-10. Unplaced scaffolds are named identically in both files, e.g. B73V4_ctg237. One other difference is that the maizegdb file includes a fasta record for the chloroplast genome (named Pt), while the fasta I used did not include plastids.
# Apologies that there are likely small differences in chromosome sequence names throughout files in this repository, like Chr1 vs chr1 vs 1!
