
### Base composition 

This counts 'methylatable' sequences in each TE and flanking region, by looking for CG, CHG, and CHH sites.

`get_base_composition_te_and_flank.sh` generates the files:
  - `B73v4.filtered_disjoined_TEs.ecology.cytpatterns.out` which is each TE
  - `B73v4.filtered_disjoined_TEs.flank.cytpatterns.out` which is each TE's 1kb flanking regions
  - `B73v4.entirechr.cytpatterns.out` which is the base composition of each chromosome (to calculate what the entire genome looks like)

`summarize_methylatable_te_flank.R` summarizes the bedtools output, writes 
- `B73_TE_methylatable.txt` with columns for `percGC`, `nCG`, `nCHG`, and `nCHH` across the TE. Slightly confusingly, the C contexts starting with `n` are the proportion of all doublets or triplets in the TE that are in the given context (e.g. )
- `B73_TE_methylatable.flank.txt` with columns for `percGC`, `nCG`, `nCHG`, and `nCHH` across the 1kb flanking region. These column names get `1kb_flank` appended in `../age_model/` to become `percGC_1kbflank`, `nCG_1kbflank`, `nCHG_1kbflank`, and `nCHH_1kbflank` 
- `B73.basecomp_genomewide.txt` for genomewide `percGC`, `nCG`, `nCHG`, and `nCHH`.

