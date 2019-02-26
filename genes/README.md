## Closest genes

`get_closest_gene.R` finds the genes closest to each TE.

*Note* that this is actually distance to an exon, so intronic TEs will be closest to an exon of the gene they are found within.
A zero distance thus means overlapping with an exon. 

Generates `B73v4_closest_gene.txt`, which includes columns for each TE copy of:

1. closest, ignoring strand - this reflects the minimum distance from this TE to a gene
  - distance (`closest`)
  - gene name (`closestgene`)
  - gene type (`closestgenetype`)

### get expression of nearby genes
- `closest_gene_expr_walley_atlas.R` uses Walley to get expression in each tissue, summarized as mean across replicates.
  - outputs `B73_closest_gene_expression.maizegdbWalley.txt`


