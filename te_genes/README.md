
identifies TE genes within each TE model, and generates gff3's with these genes, LTR, TIR, and TSD positions defined as children of the Parent TEID

- `orfs/` looks for and reports longest orfs in TE model
- `proteins/` runs hmmer using TE profile hmms to identify coding regions within TE model
- `ltr_protein_domains/` uses ltrdigest output to identify each domain of an LTR retro
