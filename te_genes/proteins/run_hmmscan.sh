#!/bin/bash -login

module load hmmer
hmmpress GenProp1044_transposon_components.hmm 
hmmscan -o test_allteorf.stdout --tblout test_allteorf.out GenProp1044_transposon_components.hmm ../orfs/.transdecoder.pep
