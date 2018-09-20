
source('../GenomeInfo.sh')


## find overlaps with TEs
bedtools intersect -names AP RP -wo -a $TEDISJOINED -b AP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed RP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed > B73v4_MNaseHS.bedout.txt

### make a file that includes just the 1kb region flanking each TE
FLANKINGFILE=


## get overlaps with flanking regions
bedtools intersect -names AP RP -wo -a $FLANKINGFILE -b AP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed RP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed > B73v4_MNaseHS.flank.bedout.txt



