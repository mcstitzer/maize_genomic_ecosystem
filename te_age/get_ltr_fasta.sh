
## these are copied from b73_ecology

zgrep "long_terminal_repeat" ../B73V4.pseudomolecule.ltrdigest.ALLsubtracted.MCSnames.20May2017.gff3.gz > B73V4.pseudomolecule.ltrdigest.ALLsubtracted.MCSnames.20May2017.ltrOnly.gff3

cat ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/ltrdigest-redo/B73V4.pseudomolecule.ltrdigest_5ltr.fas ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract/subtract*/ltrdigest/B73V4*_5ltr.fas ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/ltrdigest_unplaced/B73V4.unplaced.ltrdigest_5ltr.fas ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract_unplaced/subtract*/ltrdigest/B73*_5ltr.fas > B73V4.both_pseudo_AND_unplaced.all_5ltr.fa


## ugh, need to make sure that the pervasive ASCII issue is fixed in these three prime LTRs!
cat ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/ltrdigest-redo/B73V4.pseudomolecule.ltrdigest_3ltr.fas ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract/subtract*/ltrdigest/B73V4*_3ltr.fas ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/ltrdigest_unplaced/B73V4.unplaced.ltrdigest_3ltr.fas ~/projects/agpv4_te_annotation/ncbi_pseudomolecule/ltr/mask_subtract_unplaced/subtract*/ltrdigest/B73*_3ltr.fas > B73V4.both_pseudo_AND_unplaced.all_3ltr.fa


srun -p high --time=10-00:00 python switch_fasta_names.py B73V4.both_pseudo_AND_unplaced.all_5ltr.fa B73V4.both_pseudo_AND_unplaced.origname.newNamewFam.tab > B73V4.5ltr.newName.fa

srun -p high --time=10-00:00 python switch_fasta_names.py B73V4.both_pseudo_AND_unplaced.all_3ltr.fa B73V4.both_pseudo_AND_unplaced.origname.newNamewFam.tab > B73V4.3ltr.newName.fa
