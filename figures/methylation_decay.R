library(data.table)
library(dplyr)
library(cowplot)
library(stringr)
library(plyr)

source('../GenomeInfo.R')
source('color_palette.R')

anther=fread()
earshoot=fread()
h3k9me2=fread()


allte=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')
ind=merge(basecomp, basecomp.flank, all=T)
ind=merge(ind, diversity, all=T)
ind=merge(ind, mnase, all=T)
ind=merge(ind, allte, all=T)
ind=ind[ind$tebp>=50,] ## after disjoining, some TEs are too short to be real :( - be sure to add this to all figures!!!


anthers=colnames(a)[which(grepl('\\_anther', colnames(a)))]
SAMs=colnames(a)[which(grepl('\\_SAM', colnames(a)))]
earshoots=colnames(a)[which(grepl('\\_earshoot', colnames(a)))]
seedlingleafs=colnames(a)[which(grepl('\\_all3', colnames(a)))]
h3k9s=colnames(a)[which(grepl('h3k9', colnames(a)))]
