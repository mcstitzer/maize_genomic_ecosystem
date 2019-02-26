determines the age of each TE copy


##
- `get_fasta_of_each_TE.sh` extracts a fasta of each TE, using `collapse_chromosomes.py`
- `ltr_align/` contains code to align 5' and 3' LTRs of each copy with MAFFT, calculate pairwise distance using ape in R
- `tbl_align/` contains code to align TEs within superfamilies with MAFFT, build a tree with fasttree, and extract terminal branch length using ape in R
