
## need GENOMEFA, TEDISJOINED
TEDISJOINED=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73.structuralTEv2.2018-08-31.disjoined.gff3.gz
GENOMEFA=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73V4.both_pseudo_AND_unplaced.fa
TEFILE=/home/mstitzer/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73.structuralTEv2.2018-08-31.gff3.gz


bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CG' > B73v4.filtered_disjoined_TEs.ecology.cg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAG' > B73v4.filtered_disjoined_TEs.ecology.cag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTG' > B73v4.filtered_disjoined_TEs.ecology.ctg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCG' > B73v4.filtered_disjoined_TEs.ecology.ccg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAA' > B73v4.filtered_disjoined_TEs.ecology.caa.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAT' > B73v4.filtered_disjoined_TEs.ecology.cat.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CAC' > B73v4.filtered_disjoined_TEs.ecology.cac.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTA' > B73v4.filtered_disjoined_TEs.ecology.cta.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTT' > B73v4.filtered_disjoined_TEs.ecology.ctt.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CTC' > B73v4.filtered_disjoined_TEs.ecology.ctc.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCC' > B73v4.filtered_disjoined_TEs.ecology.ccc.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCA' > B73v4.filtered_disjoined_TEs.ecology.cca.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'CCT' > B73v4.filtered_disjoined_TEs.ecology.cct.out

### have to reverse complement everything for non-symmetric CHH sites because it is going to simply find the pattern
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'TTG' > B73v4.filtered_disjoined_TEs.ecology.ttg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'ATG' > B73v4.filtered_disjoined_TEs.ecology.atg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'GTG' > B73v4.filtered_disjoined_TEs.ecology.gtg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'TAG' > B73v4.filtered_disjoined_TEs.ecology.tag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'AAG' > B73v4.filtered_disjoined_TEs.ecology.aag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'GAG' > B73v4.filtered_disjoined_TEs.ecology.gag.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'GGG' > B73v4.filtered_disjoined_TEs.ecology.ggg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'TGG' > B73v4.filtered_disjoined_TEs.ecology.tgg.out
bedtools nuc -fi $GENOMEFA -bed $TEDISJOINED -pattern 'AGG' > B73v4.filtered_disjoined_TEs.ecology.agg.out



paste <(cat B73v4.filtered_disjoined_TEs.ecology.cg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.cag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ctg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ccg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.caa.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.cat.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.cac.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.cta.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ctt.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ctc.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ccc.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.cca.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.cct.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ttg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.atg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.gtg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.tag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.aag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.gag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.ggg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.tgg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.ecology.agg.out ) > B73v4.filtered_disjoined_TEs.ecology.cytpatterns.out


awk '{print $1"\t1\t"$2}' $GENOMEFA.fai > B73v4.entirechr.bed

bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CG' > B73v4.entirechr.cg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAG' > B73v4.entirechr.cag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTG' > B73v4.entirechr.ctg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCG' > B73v4.entirechr.ccg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAA' > B73v4.entirechr.caa.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAT' > B73v4.entirechr.cat.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CAC' > B73v4.entirechr.cac.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTA' > B73v4.entirechr.cta.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTT' > B73v4.entirechr.ctt.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CTC' > B73v4.entirechr.ctc.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCC' > B73v4.entirechr.ccc.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCA' > B73v4.entirechr.cca.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'CCT' > B73v4.entirechr.cct.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'TTG' > B73v4.entirechr.ttg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'ATG' > B73v4.entirechr.atg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'GTG' > B73v4.entirechr.gtg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'TAG' > B73v4.entirechr.tag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'AAG' > B73v4.entirechr.aag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'GAG' > B73v4.entirechr.gag.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'GGG' > B73v4.entirechr.ggg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'TGG' > B73v4.entirechr.tgg.out
bedtools nuc -fi $GENOMEFA -bed B73v4.entirechr.bed -pattern 'AGG' > B73v4.entirechr.agg.out


paste <(cat B73v4.entirechr.cg.out ) <(cut -f13 B73v4.entirechr.cag.out ) <(cut -f13 B73v4.entirechr.ctg.out ) <(cut -f13 B73v4.entirechr.ccg.out ) <(cut -f13 B73v4.entirechr.caa.out ) <(cut -f13 B73v4.entirechr.cat.out ) <(cut -f13 B73v4.entirechr.cac.out ) <(cut -f13 B73v4.entirechr.cta.out ) <(cut -f13 B73v4.entirechr.ctt.out ) <(cut -f13 B73v4.entirechr.ctc.out ) <(cut -f13 B73v4.entirechr.ccc.out ) <(cut -f13 B73v4.entirechr.cca.out ) <(cut -f13 B73v4.entirechr.cct.out ) <(cut -f13 B73v4.entirechr.ttg.out ) <(cut -f13 B73v4.entirechr.atg.out ) <(cut -f13 B73v4.entirechr.gtg.out ) <(cut -f13 B73v4.entirechr.tag.out ) <(cut -f13 B73v4.entirechr.aag.out ) <(cut -f13 B73v4.entirechr.gag.out ) <(cut -f13 B73v4.entirechr.ggg.out ) <(cut -f13 B73v4.entirechr.tgg.out ) <(cut -f13 B73v4.entirechr.agg.out ) > B73v4.entirechr.cytpatterns.out

############################################################
## FLANK
############################################################

if [ ! -f B73v4.Zm00001d.allTE.FLANK1kbEACH.gff3 ]
then
bedtools flank -b 1000 -i ../tes/B73v4.Zm00001d.allTE.gff3 -g ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/B73V4.both_pseudo_AND_unplaced.fa.fai > B73v4.Zm00001d.allTE.FLANK1kbEACH.gff3 
fi

FLANKGFF=B73v4.Zm00001d.allTE.FLANK1kbEACH.gff3

bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CG' > B73v4.filtered_disjoined_TEs.flank.cg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAG' > B73v4.filtered_disjoined_TEs.flank.cag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTG' > B73v4.filtered_disjoined_TEs.flank.ctg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCG' > B73v4.filtered_disjoined_TEs.flank.ccg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAA' > B73v4.filtered_disjoined_TEs.flank.caa.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAT' > B73v4.filtered_disjoined_TEs.flank.cat.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CAC' > B73v4.filtered_disjoined_TEs.flank.cac.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTA' > B73v4.filtered_disjoined_TEs.flank.cta.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTT' > B73v4.filtered_disjoined_TEs.flank.ctt.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CTC' > B73v4.filtered_disjoined_TEs.flank.ctc.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCC' > B73v4.filtered_disjoined_TEs.flank.ccc.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCA' > B73v4.filtered_disjoined_TEs.flank.cca.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'CCT' > B73v4.filtered_disjoined_TEs.flank.cct.out

### have to reverse complement everything for non-symmetric CHH sites because it is going to simply find the pattern
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'TTG' > B73v4.filtered_disjoined_TEs.flank.ttg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'ATG' > B73v4.filtered_disjoined_TEs.flank.atg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'GTG' > B73v4.filtered_disjoined_TEs.flank.gtg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'TAG' > B73v4.filtered_disjoined_TEs.flank.tag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'AAG' > B73v4.filtered_disjoined_TEs.flank.aag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'GAG' > B73v4.filtered_disjoined_TEs.flank.gag.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'GGG' > B73v4.filtered_disjoined_TEs.flank.ggg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'TGG' > B73v4.filtered_disjoined_TEs.flank.tgg.out
bedtools nuc -fi $GENOMEFA -bed $FLANKGFF -pattern 'AGG' > B73v4.filtered_disjoined_TEs.flank.agg.out



paste <(cat B73v4.filtered_disjoined_TEs.flank.cg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.cag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ctg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ccg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.caa.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.cat.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.cac.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.cta.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ctt.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ctc.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ccc.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.cca.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.cct.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ttg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.atg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.gtg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.tag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.aag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.gag.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.ggg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.tgg.out ) <(cut -f19 B73v4.filtered_disjoined_TEs.flank.agg.out ) > B73v4.filtered_disjoined_TEs.flank.cytpatterns.out
