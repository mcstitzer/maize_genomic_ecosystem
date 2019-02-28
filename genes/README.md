## Closest genes

`get_closest_gene.R` finds the genes closest to each TE.

*Note* that this is actually distance to an exon, so intronic TEs will be closest to an exon of the gene they are found within.
A zero distance thus means overlapping with an exon. 

Generates `B73v4_closest_gene.txt`, which includes columns for each TE copy of:

1. closest, ignoring strand - this reflects the minimum distance from this TE to a gene
  - distance (`closest`)
  - gene name (`closestgene`), AGPv4 names
  - gene type (`closestgenetype`), e.g. 5' UTR or CDS

### get expression of nearby genes
- `closest_gene_expr_walley_atlas.R` uses Walley to get expression in each tissue, summarized as mean across replicates.
  - outputs `B73_closest_gene_expression.maizegdbWalley.txt` with columns denoting the expression level of the closest gene
      - `gene_6_7_internode`    
      - `gene_7_8_internode`    
      - `gene_Vegetative_Meristem_16_19_Day`"    
      - `gene_Ear_Primordium_2_4_mm`
      - `gene_Ear_Primordium_6_8_mm`
      - `gene_Embryo_20_DAP`
      - `gene_Embryo_38_DAP`
      - `gene_Endosperm_12_DAP`
      - `gene_Endosperm_Crown_27_DAP`
      - `gene_Germinatin_Kernels_2_DAI`
      - `gene_Pericarp_Aleurone_27_DAP`
      - `gene_Leaf_Zone_1__Symmetrical_`
      - `gene_Leaf_Zone_2__Stomatal_`
      - `gene_Leaf_Zone_3__Growth_`
      - `gene_Mature_Leaf_8`
      - `gene_Primary_Root_5_Days`
      - `gene_Root___Cortex_5_Days`
      - `gene_Root___Elongation_Zone_5_Days`
      - `gene_Root___Meristem_Zone_5_Days`
      - `gene_Secondary_Root_7_8_Days`
      - `gene_B73_Mature_Pollen`
      - `gene_Female_Spikelet_Collected_on_day_as_silk`
      - `gene_Silk`
      - `gene_coefvar`
      - `gene_median`
      - `gene_tau`
      - `famsize`

Warning: *While these tissues are identical to those in ../te_expression/, their names differ slightly*
