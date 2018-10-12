#!/bin/bash -login

SUP=$1

~/software/bin/mafft --memsavetree --adjustdirection ${SUP}.fa > ${SUP}.aln.fa
~/software/FastTreeMP -nt -gtr < ${SUP}.aln.fa > ${SUP}.tre


