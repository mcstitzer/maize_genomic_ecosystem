
source ../GenomeInfo.R ## even though this is an R file, it sets up TEDISJOINED, GENOMEFA, and TEFILE in a way that bash can handle
## need GENOMEFA, TEDISJOINED
#TEDISJOINED=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73.structuralTEv2.2018-08-31.disjoined.gff3.gz
#GENOMEFA=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73V4.both_pseudo_AND_unplaced.fa
#TEFILE=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73.structuralTEv2.2018-08-31.gff3.gz


bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CG' > ${GENOME}.filtered_disjoined_TEs.ecology.cg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAG' > ${GENOME}.filtered_disjoined_TEs.ecology.cag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTG' > ${GENOME}.filtered_disjoined_TEs.ecology.ctg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCG' > ${GENOME}.filtered_disjoined_TEs.ecology.ccg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAA' > ${GENOME}.filtered_disjoined_TEs.ecology.caa.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAT' > ${GENOME}.filtered_disjoined_TEs.ecology.cat.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAC' > ${GENOME}.filtered_disjoined_TEs.ecology.cac.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTA' > ${GENOME}.filtered_disjoined_TEs.ecology.cta.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTT' > ${GENOME}.filtered_disjoined_TEs.ecology.ctt.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTC' > ${GENOME}.filtered_disjoined_TEs.ecology.ctc.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCC' > ${GENOME}.filtered_disjoined_TEs.ecology.ccc.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCA' > ${GENOME}.filtered_disjoined_TEs.ecology.cca.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCT' > ${GENOME}.filtered_disjoined_TEs.ecology.cct.out

### have to reverse complement everything for non-symmetric CHH sites because it is going to simply find the pattern
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'TTG' > ${GENOME}.filtered_disjoined_TEs.ecology.ttg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'ATG' > ${GENOME}.filtered_disjoined_TEs.ecology.atg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'GTG' > ${GENOME}.filtered_disjoined_TEs.ecology.gtg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'TAG' > ${GENOME}.filtered_disjoined_TEs.ecology.tag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'AAG' > ${GENOME}.filtered_disjoined_TEs.ecology.aag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'GAG' > ${GENOME}.filtered_disjoined_TEs.ecology.gag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'GGG' > ${GENOME}.filtered_disjoined_TEs.ecology.ggg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'TGG' > ${GENOME}.filtered_disjoined_TEs.ecology.tgg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'AGG' > ${GENOME}.filtered_disjoined_TEs.ecology.agg.out



paste <(cat B73v4.filtered_disjoined_TEs.ecology.cg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.cag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ctg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ccg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.caa.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.cat.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.cac.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.cta.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ctt.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ctc.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ccc.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.cca.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.cct.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ttg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.atg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.gtg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.tag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.aag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.gag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.ggg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.tgg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.ecology.agg.out ) > ${GENOME}.filtered_disjoined_TEs.ecology.cytpatterns.out


awk '{print $1"\t1\t"$2}' $GENOMEFA.fai > ${GENOME}.entirechr.bed

bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CG' > ${GENOME}.entirechr.cg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAG' > ${GENOME}.entirechr.cag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTG' > ${GENOME}.entirechr.ctg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCG' > ${GENOME}.entirechr.ccg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAA' > ${GENOME}.entirechr.caa.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAT' > ${GENOME}.entirechr.cat.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAC' > ${GENOME}.entirechr.cac.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTA' > ${GENOME}.entirechr.cta.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTT' > ${GENOME}.entirechr.ctt.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTC' > ${GENOME}.entirechr.ctc.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCC' > ${GENOME}.entirechr.ccc.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCA' > ${GENOME}.entirechr.cca.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCT' > ${GENOME}.entirechr.cct.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'TTG' > ${GENOME}.entirechr.ttg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'ATG' > ${GENOME}.entirechr.atg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'GTG' > ${GENOME}.entirechr.gtg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'TAG' > ${GENOME}.entirechr.tag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'AAG' > ${GENOME}.entirechr.aag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'GAG' > ${GENOME}.entirechr.gag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'GGG' > ${GENOME}.entirechr.ggg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'TGG' > ${GENOME}.entirechr.tgg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'AGG' > ${GENOME}.entirechr.agg.out


paste <(cat ${GENOME}.entirechr.cg.out ) <(cut -f13 ${GENOME}.entirechr.cag.out ) <(cut -f13 ${GENOME}.entirechr.ctg.out ) <(cut -f13 ${GENOME}.entirechr.ccg.out ) <(cut -f13 ${GENOME}.entirechr.caa.out ) <(cut -f13 ${GENOME}.entirechr.cat.out ) <(cut -f13 ${GENOME}.entirechr.cac.out ) <(cut -f13 ${GENOME}.entirechr.cta.out ) <(cut -f13 ${GENOME}.entirechr.ctt.out ) <(cut -f13 ${GENOME}.entirechr.ctc.out ) <(cut -f13 ${GENOME}.entirechr.ccc.out ) <(cut -f13 ${GENOME}.entirechr.cca.out ) <(cut -f13 ${GENOME}.entirechr.cct.out ) <(cut -f13 ${GENOME}.entirechr.ttg.out ) <(cut -f13 ${GENOME}.entirechr.atg.out ) <(cut -f13 ${GENOME}.entirechr.gtg.out ) <(cut -f13 ${GENOME}.entirechr.tag.out ) <(cut -f13 ${GENOME}.entirechr.aag.out ) <(cut -f13 ${GENOME}.entirechr.gag.out ) <(cut -f13 ${GENOME}.entirechr.ggg.out ) <(cut -f13 ${GENOME}.entirechr.tgg.out ) <(cut -f13 ${GENOME}.entirechr.agg.out ) > ${GENOME}.entirechr.cytpatterns.out

############################################################
## FLANK
############################################################

if [ ! -f ${GENOME}.allTE.FLANK1kbEACH.gff3 ]
then
bedtools flank -b 1000 -i $TEFILE -g ${GENOMEFA}.fai > ${GENOME}.allTE.FLANK1kbEACH.gff3 
fi

FLANKGFF=${GENOME}.allTE.FLANK1kbEACH.gff3

bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CG' > ${GENOME}.filtered_disjoined_TEs.flank.cg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAG' > ${GENOME}.filtered_disjoined_TEs.flank.cag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTG' > ${GENOME}.filtered_disjoined_TEs.flank.ctg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCG' > ${GENOME}.filtered_disjoined_TEs.flank.ccg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAA' > ${GENOME}.filtered_disjoined_TEs.flank.caa.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAT' > ${GENOME}.filtered_disjoined_TEs.flank.cat.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAC' > ${GENOME}.filtered_disjoined_TEs.flank.cac.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTA' > ${GENOME}.filtered_disjoined_TEs.flank.cta.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTT' > ${GENOME}.filtered_disjoined_TEs.flank.ctt.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTC' > ${GENOME}.filtered_disjoined_TEs.flank.ctc.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCC' > ${GENOME}.filtered_disjoined_TEs.flank.ccc.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCA' > ${GENOME}.filtered_disjoined_TEs.flank.cca.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCT' > ${GENOME}.filtered_disjoined_TEs.flank.cct.out

### have to reverse complement everything for non-symmetric CHH sites because it is going to simply find the pattern
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'TTG' > ${GENOME}.filtered_disjoined_TEs.flank.ttg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'ATG' > ${GENOME}.filtered_disjoined_TEs.flank.atg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'GTG' > ${GENOME}.filtered_disjoined_TEs.flank.gtg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'TAG' > ${GENOME}.filtered_disjoined_TEs.flank.tag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'AAG' > ${GENOME}.filtered_disjoined_TEs.flank.aag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'GAG' > ${GENOME}.filtered_disjoined_TEs.flank.gag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'GGG' > ${GENOME}.filtered_disjoined_TEs.flank.ggg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'TGG' > ${GENOME}.filtered_disjoined_TEs.flank.tgg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'AGG' > ${GENOME}.filtered_disjoined_TEs.flank.agg.out



paste <(cat ${GENOME}.filtered_disjoined_TEs.flank.cg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.cag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ctg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ccg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.caa.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.cat.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.cac.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.cta.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ctt.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ctc.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ccc.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.cca.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.cct.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ttg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.atg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.gtg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.tag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.aag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.gag.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.ggg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.tgg.out ) <(cut -f19 ${GENOME}.filtered_disjoined_TEs.flank.agg.out ) > ${GENOME}.filtered_disjoined_TEs.flank.cytpatterns.out

mkdir -p individual_patterns

mv *filtered_disjoined_TEs* individual_patterns
mv individual_patterns/*.cytpatterns.out .

