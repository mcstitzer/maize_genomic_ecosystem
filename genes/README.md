## Closest genes

`get_closest_gene.R` finds the genes closest to each TE, in four different ways.

These numbers are calculated from the TE's point of view, but it may also be interesting to think about the gene's point of view for other analyses.

*Note* that this is actually distance to an exon, so intronic TEs will be closest to an exon of the gene they are found within.
A zero distance thus means overlapping with an exon. 


Generates `B73v4_closest_gene.txt`, which includes columns for each TE copy of:

1. closest, ignoring strand - this reflects the minimum distance from this TE to a gene
  - distance (`closest`)
  - gene name (`closestgene`)
  - gene type (`closestgenetype`)
2. closest, same strand - this reflects the minimum distance from this TE to a gene on the same strand
  - *Note* that if a TE is unstranded (like, for example, a nonautonomous MITE we can't assign a direction), this will return the same result as 1.
  - distance (closest.samestrand)
  - gene name (closestgene.samestrand)
  - gene type (closestgenetype.samestrand)
3. closest, upstream - this reflects the nearest gene upstream of the TE, on the same strand. 
  - distance (closest.upstream)
  - gene name (closestgene.upstream)
  - gene type (closestgenetype.upstream)
4. closest, downstream - this reflects the nearest gene downstream of the TE, on the same strand. 
  - distance (closest.downstream)
  - gene name (closestgene.downstream)
  - gene type (closestgenetype.downstream)





#### not currently implemented (but I think I have code somewhere)
- head to head orientation (promoters of TE and gene are on opposite strands, transcription proceeds away from each other)
- tail to tail orientation (direction of transcription of TE crashes into direction of transcription of gene)


### get expression of nearby genes
- `closest_gene_expression.R` uses Sarah's gene expression counts to get expression in each tissue, summarized as mean across replicates.
  - outputs `B73v4_closest_gene_expression.txt`
