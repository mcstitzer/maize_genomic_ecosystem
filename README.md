# maize_genomic_ecosystem

`GenomeInfo.R` is sourced by R scripts in each directory, provides consistent name variables, and contains the full paths to genome fastas and TE annotations.


- `subgenomes/` assigns each TE to a subgenome
- `genes/` finds closest genes, and their expression across developmental atlases
- `mnase/` finds overlap of TEs and flank to MNase hypersensitive regions
- `methylation/` finds DNA and histone methylation of TEs and flanking regions
- `base_composition/` finds methylatable sequence in TE and flanking regions

## these ones need to be done

- `te_genes/` identifies TE genes within each TE model, and generates gff3's with these genes, LTR, TIR, and TSD positions defined as children of the Parent TEID
- `te_age/` determines the age of each TE copy
- `recombination/` identifies the recombinational environment the TE is found in
- `te_expression/` finds the expression level of each TE family and TE across developmental atlases
- `diversity/` finds number of segregating sites within each TE copy, and flanking sequence

