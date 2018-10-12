#!/bin/bash -login




## get the hmms for transposase
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam//releases/Pfam28.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
grep -f GenProp1044_transposon_components.acc.txt Pfam-A.hmm | awk '{print $2}' > GenProp1044_transposon_components.acc.Pfam-A.txt
hmmfetch -f Pfam-A.hmm GenProp1044_transposon_components.acc.Pfam-A.txt > GenProp1044_transposon_components.hmm
rm Pfam-A.hmm
