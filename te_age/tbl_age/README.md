#### calculate terminal branch length ages

1. get the first or last 1000 bp for the superfamilies that have very long TEs and so many TEs it's not computationally feasible to align each copy
  - `get_terminal_bp_sups.sh`
  - this gets terminal fasta for RLC, RLG, RLX, and DHH

1. submit each superfamily or fasta prefix to be aligned and have tree built
  - e.g. for entire TE length of a superfamily `align_TE_families.sh DTA`
  - e.g. for 1000bp of a TE superfamily `align_TE_families.sh RLC.1000`
  - each outputs a .tre and .aln.fa file for the superfamily
  
3. get terminal branch lengths for each copy
  - `Rscript summarize_trees.R` 
