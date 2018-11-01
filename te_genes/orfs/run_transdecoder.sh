#!/bin/bash -login

~/software/TransDecoder-3.0.1/TransDecoder.LongOrfs -m 20 -t ../../te_age/B73.eachTE.fa

~/software/TransDecoder-3.0.1/TransDecoder.Predict -t ../../te_age/B73.eachTE.fa --single_best_orf

Rscript filter_transdecoder.R
