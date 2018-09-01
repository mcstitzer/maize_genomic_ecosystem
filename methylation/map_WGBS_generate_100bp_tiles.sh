
module load picardtools/1.107 ### this is an earlier version where individual jar files still exist (not all merged into picardtools.jar)
module load bamtools

CPU=10

##### manually looped through each of these
#TISSUE=B73_all3
#FQ1=B73_all3_R1_val_1.fq
#FQ2=B73_all3_R2_val_2.fq
#TISSUE=B73_anther
FQ1=B73_anther_WGBS_R1_val_1.fq
FQ2=B73_anther_WGBS_R2_val_2.fq
TISSUE=B73_earshoot
#FQ1=B73_earshoot_WGBS_R1_val_1.fq
#FQ2=B73_earshoot_WGBS_R2_val_2.fq
#TISSUE=B73_flag_leaf
#FQ1=B73_flag_leaf_WGBS_R1_val_1.fq
#FQ2=B73_flag_leaf_WGBS_R2_val_2.fq
#TISSUE=B73_SAM
#FQ1=B73_SAM_WGBS_R1_val_1.fq
#FQ2=B73_SAM_WGBS_R2_val_2.fq


GENOME=B73V4.both_pseudo_AND_unplaced.chrprefix.fa

mkdir -p $TISSUE

# 1, align adapter trimmed datasets to B73 genome
bsmap -a $FQ1 -b $FQ2 -d $GENOME -o $TISSUE/${TISSUE}.wgbs.bam -v 5 -r 0 -p $CPU -q 20 -A AGATCGGAAGAGCGGTTCAGCAGGAATGCCG

# 2a, reomve PCR duplicates, must be sorted by coordinate
java -jar /share/apps/picardtools-1.107/SortSam.jar INPUT=$TISSUE/${TISSUE}.wgbs.bam OUTPUT=$TISSUE/${TISSUE}.wgbs.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT ## adding validation stringency because contigs are having problems with boundaries, i think reads mapping near an edge are considered on different chr/contigs? 
java -jar /share/apps/picardtools-1.107/MarkDuplicates.jar I=$TISSUE/${TISSUE}.wgbs.sorted.bam O=$TISSUE/${TISSUE}.wgbs.sorted.MarkDup.bam METRICS_FILE=$TISSUE/${TISSUE}.wgbs.sorted.MarkDup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT ### need to add it here too.

# 2b, keep properly paired reads
bamtools filter -isProperPair true -in $TISSUE/${TISSUE}.wgbs.sorted.MarkDup.bam -out $TISSUE/${TISSUE}.wgbs.sorted.MarkDup.pairs.bam

# 2c, clip overlapping reads
bam clipOverlap --in $TISSUE/${TISSUE}.wgbs.sorted.MarkDup.pairs.bam --out $TISSUE/${TISSUE}.wgbs.sorted.MarkDup.pairs.clipOverlap.bam

# 3, extract DNA methylation at each cysotine
python bsmap-2.74/methratio.py -o $TISSUE/${TISSUE}.wgbs.clipOverlap.bsmap274 -d $GENOME -u -z -r $TISSUE/${TISSUE}.wgbs.sorted.MarkDup.pairs.clipOverlap.bam

# 4a, get 100-bp tiles
for file in $TISSUE/${TISSUE}.wgbs.clipOverlap.bsmap274
do
	awk '{if(($3=="-" && $4~/^.CG../ ) || ($3=="+" &&  $4~/^..CG./)) print $1"\t"$2"\t"$3"\t""CG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^C[AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT]G/)) print $1"\t"$2"\t"$3"\t""CHG""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else if(($3=="-" && $4~/^[AGT][AGT]G../ ) || ($3=="+" &&  $4~/^..C[ACT][ACT]/)) print $1"\t"$2"\t"$3"\t""CHH""\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12; else print $0}' $file > $file"_BSMAP_out.txt"
	perl met_context_window.pl $file"_BSMAP_out.txt" 100           
done

# 4b, merge files based on coordinates

### generate window file
samtools faidx $GENOME
bedtools makewindows -g ${GENOME}.fai  -w 100 > ZmV4_w100.bed

### summarize data
for i in */*.wgbs.clipOverlap.bsmap274
do
	echo $i
	for j in CG CHG CHH
	do
		awk -F "\t" 'NR==FNR {f1[$1"\t"$2"\t"$3] = $4"\t"$5"\t"$6; next}  {idx=$1"\t"$2"\t"$3; if(idx in f1) print $0"\t"f1[idx]; else print $0"\t""NA""\t""NA""\t""NA"}' $i"_BSMAP_out.txt.100."$j".bed" ZmV4_w100.bed > temp.txt 
		mv temp.txt ZmV4_w100.bed
done
done

### change name, add header
mv ZmV4_w100.bed ZmV4_w100.fivetissues.bed
header="chr\tstart\tend"
for i in */*wgbs.clipOverlap.bsmap274    ## ex. B73_all3/B73_all3.wgbs.clipOverlap.bsmap274
do
	for j in CG CHG CHH
	do
		tissue=$( basename $i ".wgbs.clipOverlap.bsmap274")
		header="${header}\tCT_${tissue}_${j}"
		header="${header}\tC_${tissue}_${j}"
		header="${header}\tCratio_${tissue}_${j}"
	done
done

printf "${header}\n" > tempheader.bed
cat ZmV4_w100.fivetissues.bed >> tempheader.bed
mv tempheader.bed ZmV4_w100.fivetissues.bed

# file ZmV2_w100.bed shows non-overlapping sliding 100-bp tile of maize genome, and has the following format 
#	chr	start	end
#	chr1	0	100
#	chr1	101	200
