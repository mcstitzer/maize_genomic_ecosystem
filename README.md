# maize_genomic_ecosystem

Scripts used in ["The Genomic Ecosystem of Transposable Elements in Maize," Stitzer et al., 2021](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009768)

Interactive distributions of raw data of each TE family at https://mcstitzer.shinyapps.io/maize_te_families/

`GenomeInfo.R` is sourced by R scripts in each directory, provides consistent name variables, and contains the full paths to genome fastas and TE annotations (most of the time - some paths are hard coded in the following directories)

- `age_model/` predicts age using a random forest model
- `base_composition/` finds methylatable sequence in TE and flanking regions
- `diversity/` finds number of segregating sites within each TE copy, and flanking sequence
- `figures/` has code for generating manuscript figures
- `genes/` finds closest genes, and their expression across developmental atlases
- `methylation/` finds DNA and histone methylation of TEs and flanking regions
- `mnase/` finds overlap of TEs and flank to MNase hypersensitive regions
- `recombination/` identifies the recombinational environment the TE is found in
- `subgenomes/` assigns each TE to a subgenome
- `te_age/` determines the age of each TE copy
- `te_characteristics/` identifies features specific to the TE copy encoded within the gff, like TE length
- `te_expression/` calculates the per-copy TE expression level across developmental atlases
- `te_genes/` identifies TE genes within each TE model

