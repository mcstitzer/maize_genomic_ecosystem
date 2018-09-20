
source('../GenomeInfo.R') ## works as bash as well as R



## find overlaps with TEs
bedtools intersect -names AP RP -wo -a $TEDISJOINED -b AP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed RP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed > ${GENOME}_MNaseHS.bedout.txt

### make a file that includes just the 1kb region flanking each TE
if [ ! -f ${GENOMEFA}.fai ]
then
samtools faidx ${GENOMEFA}
fi

if [ ! -f ${GENOME}.allTE.FLANK1kbEACH.gff3 ]
then
bedtools flank -b 1000 -i $TEFILE -g ${GENOMEFA}.fai > ${GENOME}.allTE.FLANK1kbEACH.gff3 
fi

FLANKGFF=${GENOME}.allTE.FLANK1kbEACH.gff3


## get overlaps with flanking regions
bedtools intersect -names AP RP -wo -a $FLANKGFF -b AP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed RP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed > ${GENOME}_MNaseHS.flank.bedout.txt



