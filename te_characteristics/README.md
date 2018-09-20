## general info on each TE copy

`count_bp_disruptors.R` generates a file named like `B73_TE_individual_copies.2018-09-19.txt`:
  - Col10, `tebp`: the number of bp this TE contributes to the genome
  - Col11, `tespan`: the physical distance this TE measures along the genome (including bp actually contributed by TEs nested into it)
  - Col12, `pieces`: the number of pieces the TE is cut into by being disrupted by other TEs
  - Col13, `disruptor`: the number of other TEs this TE is inserted into
