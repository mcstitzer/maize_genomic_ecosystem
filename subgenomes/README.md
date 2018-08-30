## Subgenomes

Two ancestral subgenomes can be identified in maize, resulting from an allopolyploidy event following hybridization of two divergent African grasses [McCain et al., 2018](https://www.biorxiv.org/content/early/2018/06/20/352351). 

Using subgenome calls from [Jiao et al., 2017](https://www.nature.com/articles/nature22971), we find which subgenome each TE can be assigned to in `assign_TEs_to_subgenomes.R`.

- This outputs four files, with filenames like this (but including the date they were generated before the .txt):
  - `B73v4_TE_subgenome.txt` counting the subgenome for each TE
  - `B73v4_TE_subgenome.fam.txt` counting the number of copies of a family in each subgenome
  - `B73v4_TE_subgenome.sup.txt` counting the number of copies of a superfamily in each subgenome
  - `genome-wide-proportion-subgenomeA.txt` where the proportional number of bp assigned to the A subgenome is here
