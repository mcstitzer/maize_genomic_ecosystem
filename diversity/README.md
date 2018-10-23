finds number of segregating sites within each TE copy, and flanking sequence

## uses HapMap3.2.1 unimputed SNPs, ported to AGPv4

These are ten chromosome files from panzea, details [here](http://cbsusrv04.tc.cornell.edu/users/panzea/filegateway.aspx?category=Genotypes), and available here: `/iplant/home/shared/panzea/hapmap3/hmp321/unimputed/uplifted_APGv4/hmp321_agpv4_chr${i}.vcf.gz`

- `calc_seg_sites.sh` generates intermediate files that download vcf's, concatenate to bed, and overlap with TE and flank.
  - this is optimized for the UCD cluster, as it uses `/scratch/` to save disk I/O

- `summarize_seg_sites.R` generates file `GENOME.segregatingsites.TEandFlank.txt`, with columns:
  - `TEID`, `fam`, and `sup`
  - `te_bp` the bp of the TE itself
  - `segsites` the count of segregating sites in the TE
  - `segsites.bp` the proportion of sites in the TE that are segregating
  - `flank_bp` the number of flanking bases assayed (2000 bp, unless at the end of a chromosome or contig)
  - `flank_segsites` the count of segregating sites in the flanking sequence
  - `flank_segsites.bp` the proportion of sites in the flanking sequence that are segregating
