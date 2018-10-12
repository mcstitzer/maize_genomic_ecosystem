#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/b73_ecology/RawData/tes/age_recalculation
#SBATCH -o /home/mstitzer/projects/b73_ecology/slurm-log/align_consensus-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/b73_ecology/slurm-log/align_consensus-stderr-%j.txt
#SBATCH -J align_cons
set -e
set -u

module load emboss
module load mafft

mkdir -p split_fasta_5p
python split_fasta.py B73V4.5ltr.newName.fa split_fasta_5p

mkdir -p split_fasta_3p
python split_fasta.py B73V4.3ltr.newName.fa split_fasta_3p

TEID=($(awk '{print $2}' B73V4.both_pseudo_AND_unplaced.origname.newNamewFam.tab))

mkdir -p aligned_ltrs
mkdir -p consensus_ltrs

for TE in "${TEID[@]}"
do
 if [ ! -f aligned_ltrs/${TE}.aln.fa ]
 then
 if [ -f split_fasta_5p/${TE}.fa ]
  then
  cat split_fasta_5p/${TE}.fa split_fasta_3p/${TE}.fa > temp.fa
  mafft temp.fa > aligned_ltrs/${TE}.aln.fa 2> /dev/null
  consambig -name ${TE} -sequence aligned_ltrs/${TE}.aln.fa -outseq consensus_ltrs/${TE}.cons.fa 2> /dev/null
 fi
 fi
done

rm temp.fa
