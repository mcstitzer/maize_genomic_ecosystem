library(dplyr)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(plyr)

source('color_palette.R')

ind=fread('../age_model/B73.LTRAGE.allDescriptors.2019-10-22.txt')
ind$sup=as.factor(ind$sup)
ind$mya=ind$age/2/3.3e-8/1e6

imp=fread('../age_model/importance_each_variable.2019-10-23.txt')
impsort=fread('../age_model/importance_combo_variable.2019-10-23.txt')

largest10=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),10)))))

#largest10=largest10[-which(substr(names(largest10),1,3)=='RIL')[-c(1,2)]]                 
largest10
ind$largest10=ind$fam %in% names(largest10)                               
largest10=largest10[c(1:70,83:112,71:82,113:122)] ## super hard coded to get the nonLTR together!                     
largest5=unlist(unique(sapply(unique(ind$sup), function(x) rev(tail(sort(table(ind$fam[ind$sup==x & !duplicated(ind$TEID)])),5)))))
largest5
largest5=largest5[c(1:35,43:57,36:42,58:62)] ## super hard coded to get the nonLTR together!                     
ind$largest5=ind$fam %in% names(largest5)                               


## get rmse in units of mya
imp$rmseMya
imp$rmseMya=sqrt(imp$X.IncMSE/100)/2/3.3e-8/1e6
impsort$rmseMya=sqrt(impsort$sum/100)/2/3.3e-8/1e6 ## should i be doing this before summing though?

## make a color scale for importance plots!
myColors <- brewer.pal(12,"Paired")
myColors=myColors[c(1:6,8,10,11,12)]
#names(myColors) <- c('flank_base_composition', 'TE_base_composition', 'flank_methylation_mnase', 'TE_methylation_mnase',
#                     'flank_closest_gene_expression', 'TE_expression', 'TE_genes', 'TE_taxonomy', 'flank_selection', 'TE_features')
names(myColors) <- c('Flank base composition', 'TE base composition', 'Flank methylation and MNase', 'TE methylation and MNase',
                     'Flanking gene expression', 'TE expression', 'TE encoded genes', 'TE taxonomy','Flank selection', 'TE features')

colScale <- scale_colour_manual(name = "feature",values = myColors)

## and assign to categories
## assign to categories
categories=data.frame(feature=imp$feat)
categories$category=NA
categories$category[categories$feature %in% c('fam', 'sup')]='TE_taxonomy'
categories$category[categories$feature %in% c('famsize', 'pieces', 'chr', 'helLCV3p', 'helLCV5p', 'TIRlen', 'avg_ltr_length','te_bp', 'tebp', 'tespan' )]='TE_features'
categories$category[grepl('_avg_c', categories$feature) | categories$feature %in% c('shoot_prop', 'shoot_bp','root_bp', 'root_prop', 'n_shoot_hs', 'n_root_hs') ]='TE_methylation_mnase'
categories$category[categories$feature %in% c('nCG', 'nCHG', 'nCHH', 'nTG', 'percGC') ]='TE_base_composition'
categories$category[categories$feature %in% c('GAG', 'AP', 'RT', 'RNaseH', 'INT', 'ENV', 'CHR', 'auton', 'pol', 'orfAA', 'helprot', 'tirprot', 'rveprot', 'tpaseprot', 'GAGfam', 'APfam', 'INTfam', 'RTfam', 'RNaseHfam', 'ENVfam', 'CHRfam', 'polfam', 'autonfam', 'helprotfam', 'rveprotfam', 'tpaseprotfam')]='TE_genes'
categories$category[grepl('_[[:digit:]]+', categories$feature) | grepl('h3k9me2_[[:digit:]]+', categories$feature) | categories$feature %in% c('flank_h3k9me2', 'flank_cg', 'flank_chg', 'flank_chh')  | (grepl('_flank', categories$feature) &  grepl('oot', categories$feature))]='flank_methylation_mnase'
categories$category[grepl('^avg_c', categories$feature)]='flank_base_composition'
categories$category[grepl('cm', categories$feature) | grepl('segsites', categories$feature) | categories$feature %in% c('closest', 'closestgenesyntenic', 'closest.syntenic', 'ingene', 'subgenome', 'disruptor')]='flank_selection'
categories$category[grepl('^gene_', categories$feature) | grepl('^syntenicgene_', categories$feature) ]='flank_closest_gene_expression'
categories$category[grepl('^TEfam', categories$feature) | grepl('unique', categories$feature)]='TE_expression'

## new ways to get at this - quick fix :(
categories$category[categories$feature %in% c('percGC_1kbflank','nCG_1kbflank','nCHG_1kbflank','nCHH_1kbflank','flank_bp', 'nTG_1kb_flank')]='flank_base_composition'
categories$category[categories$feature %in% c('flank_n_root_hs','flank_root_bp','flank_root_prop','flank_n_shoot_hs','flank_shoot_bp','flank_shoot_prop')]='flank_methylation_mnase'
categories$category[categories$feature %in% c('segsites.bp', 'disruptor', 'disruptor.samefam', 'disruptor.samesup')]='TE_features'
categories[is.na(categories$category),]

## make these more readable for the figure
nicecategorynames=c('TE taxonomy', 'Flank selection', 'TE features', 'TE encoded genes', 'TE base composition', 
                    'Flank base composition', 'Flank methylation and MNase', 'TE methylation and MNase', 'Flanking gene expression',
                    'TE expression')
names(nicecategorynames)=c("TE_taxonomy", "flank_selection", "TE_features", "TE_genes", 
                          "TE_base_composition", "flank_base_composition", "flank_methylation_mnase", 
                          "TE_methylation_mnase", "flank_closest_gene_expression", "TE_expression")
categories$nicecategory=str_replace_all(categories$category, nicecategorynames)
nicefeaturenames=c('TE family', 'TE superfamily', 'Distance to closest gene', 'Distance to closest syntenic gene', 'TE length (bp)', 'TE span (bp)',
                   'Pieces TE disrupted into', 'TE disrupts another TE', 'TE disrupts a family member', 'TE disrupts a superfamily member',
                   'Is the closest gene syntenic?', 'Is the TE in a gene?', 'TE contains Helitron protein domain', 'TE contains Reverse Transcriptase protein domain', 'TE contains Transposase protein domain',
                   'TE sequence matches GAG protein', 'TE sequence matches Aspartic Proteinase', 'TE sequence matches Integrase', 'TE sequence matches Reverse Transcriptase', 'TE sequence matches RNaseH', 'TE sequence matches Envelope protein', 'TE sequence matches Chromodomain protein', 'TE contains matches to entire Polyprotein', 'TE contains matches to both GAG and Pol',
                   'Family member contains GAG', 'Family member contains Aspartic Proteinase', 'Family member contains Integrase', 'Family member contains Reverse Transcriptase', 'Family member contains RNaseH', 'Family member contains Envelope', 
                   'Family member contains Chromodomain', 'Family member contains entire Polyprotein', 'Family member contains both GAG and Pol', 'Family member contains Helitron protein', 'Family member contains Reverse Transcriptase', 'Family member contains Transposase',
                   'Length of longest ORF in AAs', 'Percent GC', 'Proportion CG methylatable', 'Proportion CHG methylatable', 'Proportion CHH methylatable', 'Proportion TG or CA (deaminated)', 'Proportion GC in 1kb flanking TE',
                   'Proportion CG methylatable in 1kb flank', 'Proportion CHG methylatable in 1kb flank', 'Proportion CHH methylatable in 1kb flank', 'Proportion TG or CA (deamination) in 1kb flank',
                   'Number root MNase hypersensitive peaks in TE', 'Root MNase hypersensitive bp in TE', 'Proportion of TE bp that are root MNase hypersensitive', 'Number shoot MNase hypersensitive peaks in TE', 
                   'Shoot MNase hypersensitive bp in TE', ## this one belongs on previous line
                   'Proportion of TE bp that are shoot MNase hypersensitive', 'Number root MNase hypersensitive peaks in 1kb flank', 'Root MNase hypersensitive bp in 1kb flank', 'Proportion of 1kb flank that are root MNase hypersensitive',
                   'Number shoot MNase hypersensitive peaks in 1kb flank', 'Shoot MNase hypersensitive bp in 1kb flank', 'Proportion of 1kb flank that are shoot MNase hypersensitive', 'Number segregating sites per bp across HapMap3',
                   'Number segregating sites per bp in 1kb flank in HapMap3', 'Average TE CG methylation in anther', 'Average TE CHG methylation in anther', 'Average TE CHH methylation in anther',
                   paste0('Average CG methylation ', seq(100,2000, by=100), ' bp away from TE in anther'),
                   paste0('Average CHG methylation ', seq(100,2000, by=100), ' bp away from TE in anther'),
                   paste0('Average CHH methylation ', seq(100,2000, by=100), ' bp away from TE in anther'),
                   'Average TE CG methylation in earshoot', 'Average TE CHG methylation in earshoot', 'Average TE CHH methylation in earshoot',
                   paste0('Average CG methylation ', seq(100,2000, by=100), ' bp away from TE in earshoot'),
                   paste0('Average CHG methylation ', seq(100,2000, by=100), ' bp away from TE in earshoot'),
                   paste0('Average CHH methylation ', seq(100,2000, by=100), ' bp away from TE in earshoot'),
                   'Average TE CG methylation in flag leaf', 'Average TE CHG methylation in flag leaf', 'Average TE CHH methylation in flag leaf',
                   paste0('Average CG methylation ', seq(100,2000, by=100), ' bp away from TE in flag leaf'),
                   paste0('Average CHG methylation ', seq(100,2000, by=100), ' bp away from TE in flag leaf'),
                   paste0('Average CHH methylation ', seq(100,2000, by=100), ' bp away from TE in flag leaf'),
                   'Average TE CG methylation in seedling', 'Average TE CHG methylation in seedling', 'Average TE CHH methylation in seedling',
                   paste0('Average CG methylation ', seq(100,2000, by=100), ' bp away from TE in seedling'),
                   paste0('Average CHG methylation ', seq(100,2000, by=100), ' bp away from TE in seedling'),
                   paste0('Average CHH methylation ', seq(100,2000, by=100), ' bp away from TE in seedling'),
                   'Average TE CG methylation in shoot apical meristem', 'Average TE CHG methylation in shoot apical meristem', 'Average TE CHH methylation in shoot apical meristem',
                   paste0('Average CG methylation ', seq(100,2000, by=100), ' bp away from TE in shoot apical meristem'),
                   paste0('Average CHG methylation ', seq(100,2000, by=100), ' bp away from TE in shoot apical meristem'),
                   paste0('Average CHH methylation ', seq(100,2000, by=100), ' bp away from TE in shoot apical meristem'),
                   'Expression of closest gene in 6th-7th internode', 'Expression of closest gene in 7th-8th internode',
                   'Expression of closest gene in vegetative meristem (16-19 days)', 'Expression of closest gene in 2-4mm ear primordium', 
                   'Expression of closest gene in 6-8mm ear primordium', 'Expression of closest gene in 20 DAP embryo', 'Expression of closest gene in 38 DAP embryo',
                   'Expression of closest gene in 12 DAP endosperm', 'Expression of closest gene in 27 DAP endosperm crown', 'Expression of closest gene in germinating kernels, 2 DAI',
                   'Expression of closest gene in pericarp/aleurone 27 DAP', 'Expression of closest gene in leaf 8 (leaf zone 1, symmetrical)', 
                   'Expression of closest gene in leaf 8 (leaf zone 2, stomatal)', 'Expression of closest gene in leaf 8 (leaf zone 3, growth)', 'Expression of closest gene in mature leaf 8',
                   'Expression of closest gene in 5 day old primary root', 'Expression of closest gene in 5 day old root cortex', 'Expression of closest gene in 5 day old root elongation zone',
                   'Expression of closest gene in 5 day old meristem zone root', 'Expression of closest gene in 7-8 day secondary root', 
                   'Expression of closest gene in mature pollen', 'Expression of closest gene in female spikelet', 
                   'Expression of closest gene in silk', 'Coefficient of variation of closest gene across tissues', 'Median expression of closest gene across tissues', 'Tau of closest gene', 
                   'Expression of closest syntenic gene in 6th-7th internode', 'Expression of closest syntenic gene in 7th-8th internode',
                   'Expression of closest syntenic gene in vegetative meristem (16-19 days)', 'Expression of closest syntenic gene in 2-4mm ear primordium', 
                   'Expression of closest syntenic gene in 6-8mm ear primordium', 'Expression of closest syntenic gene in 20 DAP embryo', 'Expression of closest syntenic gene in 38 DAP embryo',
                   'Expression of closest syntenic gene in 12 DAP endosperm', 'Expression of closest syntenic gene in 27 DAP endosperm crown', 'Expression of closest syntenic gene in germinating kernels, 2 DAI',
                   'Expression of closest syntenic gene in pericarp/aleurone 27 DAP', 'Expression of closest syntenic gene in leaf 8 (leaf zone 1, symmetrical)', 
                   'Expression of closest syntenic gene in leaf 8 (leaf zone 2, stomatal)', 'Expression of closest syntenic gene in leaf 8 (leaf zone 3, growth)', 'Expression of closest syntenic gene in mature leaf 8',
                   'Expression of closest syntenic gene in 5 day old primary root', 'Expression of closest syntenic gene in 5 day old root cortex', 'Expression of closest syntenic gene in 5 day old root elongation zone',
                   'Expression of closest syntenic gene in 5 day old meristem zone root', 'Expression of closest syntenic gene in 7-8 day secondary root', 
                   'Expression of closest syntenic gene in mature pollen', 'Expression of closest syntenic gene in female spikelet', 
                   'Expression of closest syntenic gene in silk', 'Coefficient of variation of closest syntenic gene across tissues', 'Median expression of closest syntenic gene across tissues', 'Tau of closest syntenic gene',
                   'Expression of TE family in 2-4mm ear primordium', 'Expression of TE family in 6-8mm ear primordium',
                   'Expression of TE family in 20 DAP embryo', 'Expression of TE family in 38 DAP embryo', 'Expression of TE family in 12 DAP endosperm',
                   'Expression of TE family in 27 DAP endosperm crown', 'Expression of TE family in 2 DAI germinating kernels', 'Expression of TE family in 6th-7th internode',
                   'Expression of TE family in 7th-8th internode', 'Expression of TE family in mature leaf 8', 'Expression of TE family in leaf 8 (leaf zone 3, growth)',
                   'Expression of TE family in leaf 8 (leaf zone 2, stomatal)', 'Expression of TE family in leaf 8 (leaf zone 1, symmetrical)', 'Expression of TE family in vegetative meristem (16-19 days)',
                   'Expression of TE family in mature pollen', 'Expression of TE family in 5 day old root cortex', 'Expression of TE family in 5 day old root elongation zone',
                   'Expression of TE family in 5 day old root meristem zone', 'Expression of 5 day old primary root', 'Expression of 5 day old secondary root',
                   'Expression of TE family in pericarp/aleurone 27 DAP', 'Expression of TE family in silk', 'Expression of TE family in female spikelet',
                   'Median expression of TE family across tissues', 'Median expression of each TE copy across tissues', 'Median expression per bp of each TE', 'Tau of expression of TE family',
                   'Recombination rate (cM/Mb)', 'Subgenome', 'Family size')

                   
names(nicefeaturenames)=c("fam", "sup", "closest", "closest.syntenic", "tebp", "tespan", 
"pieces", "disruptor", "disruptor.samefam", "disruptor.samesup", 
"closestgenesyntenic", "ingene", "helprot", "rveprot", "tpaseprot", 
"GAG", "AP", "INT", "RT", "RNaseH", "ENV", "CHR", "pol", "auton", 
"GAGfam", "APfam", "INTfam", "RTfam", "RNaseHfam", "ENVfam", 
"CHRfam", "polfam", "autonfam", "helprotfam", "rveprotfam", "tpaseprotfam", 
"orfAA", "percGC", "nCG", "nCHG", "nCHH", "nTG", "percGC_1kbflank", 
"nCG_1kbflank", "nCHG_1kbflank", "nCHH_1kbflank", "nTG_1kbflank", 
"n_root_hs", "root_bp", "root_prop", "n_shoot_hs", "shoot_bp", 
"shoot_prop", "flank_n_root_hs", "flank_root_bp", "flank_root_prop", 
"flank_n_shoot_hs", "flank_shoot_bp", "flank_shoot_prop", "segsites.bp", 
"flank_segsites.bp", "anther_avg_cg", "anther_avg_chg", "anther_avg_chh", 
"anther_flank_cg_100", "anther_flank_cg_200", "anther_flank_cg_300", 
"anther_flank_cg_400", "anther_flank_cg_500", "anther_flank_cg_600", 
"anther_flank_cg_700", "anther_flank_cg_800", "anther_flank_cg_900", 
"anther_flank_cg_1000", "anther_flank_cg_1100", "anther_flank_cg_1200", 
"anther_flank_cg_1300", "anther_flank_cg_1400", "anther_flank_cg_1500", 
"anther_flank_cg_1600", "anther_flank_cg_1700", "anther_flank_cg_1800", 
"anther_flank_cg_1900", "anther_flank_cg_2000", "anther_flank_chg_100", 
"anther_flank_chg_200", "anther_flank_chg_300", "anther_flank_chg_400", 
"anther_flank_chg_500", "anther_flank_chg_600", "anther_flank_chg_700", 
"anther_flank_chg_800", "anther_flank_chg_900", "anther_flank_chg_1000", 
"anther_flank_chg_1100", "anther_flank_chg_1200", "anther_flank_chg_1300", 
"anther_flank_chg_1400", "anther_flank_chg_1500", "anther_flank_chg_1600", 
"anther_flank_chg_1700", "anther_flank_chg_1800", "anther_flank_chg_1900", 
"anther_flank_chg_2000", "anther_flank_chh_100", "anther_flank_chh_200", 
"anther_flank_chh_300", "anther_flank_chh_400", "anther_flank_chh_500", 
"anther_flank_chh_600", "anther_flank_chh_700", "anther_flank_chh_800", 
"anther_flank_chh_900", "anther_flank_chh_1000", "anther_flank_chh_1100", 
"anther_flank_chh_1200", "anther_flank_chh_1300", "anther_flank_chh_1400", 
"anther_flank_chh_1500", "anther_flank_chh_1600", "anther_flank_chh_1700", 
"anther_flank_chh_1800", "anther_flank_chh_1900", "anther_flank_chh_2000", 
"earshoot_avg_cg", "earshoot_avg_chg", "earshoot_avg_chh", "earshoot_flank_cg_100", 
"earshoot_flank_cg_200", "earshoot_flank_cg_300", "earshoot_flank_cg_400", 
"earshoot_flank_cg_500", "earshoot_flank_cg_600", "earshoot_flank_cg_700", 
"earshoot_flank_cg_800", "earshoot_flank_cg_900", "earshoot_flank_cg_1000", 
"earshoot_flank_cg_1100", "earshoot_flank_cg_1200", "earshoot_flank_cg_1300", 
"earshoot_flank_cg_1400", "earshoot_flank_cg_1500", "earshoot_flank_cg_1600", 
"earshoot_flank_cg_1700", "earshoot_flank_cg_1800", "earshoot_flank_cg_1900", 
"earshoot_flank_cg_2000", "earshoot_flank_chg_100", "earshoot_flank_chg_200", 
"earshoot_flank_chg_300", "earshoot_flank_chg_400", "earshoot_flank_chg_500", 
"earshoot_flank_chg_600", "earshoot_flank_chg_700", "earshoot_flank_chg_800", 
"earshoot_flank_chg_900", "earshoot_flank_chg_1000", "earshoot_flank_chg_1100", 
"earshoot_flank_chg_1200", "earshoot_flank_chg_1300", "earshoot_flank_chg_1400", 
"earshoot_flank_chg_1500", "earshoot_flank_chg_1600", "earshoot_flank_chg_1700", 
"earshoot_flank_chg_1800", "earshoot_flank_chg_1900", "earshoot_flank_chg_2000", 
"earshoot_flank_chh_100", "earshoot_flank_chh_200", "earshoot_flank_chh_300", 
"earshoot_flank_chh_400", "earshoot_flank_chh_500", "earshoot_flank_chh_600", 
"earshoot_flank_chh_700", "earshoot_flank_chh_800", "earshoot_flank_chh_900", 
"earshoot_flank_chh_1000", "earshoot_flank_chh_1100", "earshoot_flank_chh_1200", 
"earshoot_flank_chh_1300", "earshoot_flank_chh_1400", "earshoot_flank_chh_1500", 
"earshoot_flank_chh_1600", "earshoot_flank_chh_1700", "earshoot_flank_chh_1800", 
"earshoot_flank_chh_1900", "earshoot_flank_chh_2000", "flagleaf_avg_cg", 
"flagleaf_avg_chg", "flagleaf_avg_chh", "flagleaf_flank_cg_100", 
"flagleaf_flank_cg_200", "flagleaf_flank_cg_300", "flagleaf_flank_cg_400", 
"flagleaf_flank_cg_500", "flagleaf_flank_cg_600", "flagleaf_flank_cg_700", 
"flagleaf_flank_cg_800", "flagleaf_flank_cg_900", "flagleaf_flank_cg_1000", 
"flagleaf_flank_cg_1100", "flagleaf_flank_cg_1200", "flagleaf_flank_cg_1300", 
"flagleaf_flank_cg_1400", "flagleaf_flank_cg_1500", "flagleaf_flank_cg_1600", 
"flagleaf_flank_cg_1700", "flagleaf_flank_cg_1800", "flagleaf_flank_cg_1900", 
"flagleaf_flank_cg_2000", "flagleaf_flank_chg_100", "flagleaf_flank_chg_200", 
"flagleaf_flank_chg_300", "flagleaf_flank_chg_400", "flagleaf_flank_chg_500", 
"flagleaf_flank_chg_600", "flagleaf_flank_chg_700", "flagleaf_flank_chg_800", 
"flagleaf_flank_chg_900", "flagleaf_flank_chg_1000", "flagleaf_flank_chg_1100", 
"flagleaf_flank_chg_1200", "flagleaf_flank_chg_1300", "flagleaf_flank_chg_1400", 
"flagleaf_flank_chg_1500", "flagleaf_flank_chg_1600", "flagleaf_flank_chg_1700", 
"flagleaf_flank_chg_1800", "flagleaf_flank_chg_1900", "flagleaf_flank_chg_2000", 
"flagleaf_flank_chh_100", "flagleaf_flank_chh_200", "flagleaf_flank_chh_300", 
"flagleaf_flank_chh_400", "flagleaf_flank_chh_500", "flagleaf_flank_chh_600", 
"flagleaf_flank_chh_700", "flagleaf_flank_chh_800", "flagleaf_flank_chh_900", 
"flagleaf_flank_chh_1000", "flagleaf_flank_chh_1100", "flagleaf_flank_chh_1200", 
"flagleaf_flank_chh_1300", "flagleaf_flank_chh_1400", "flagleaf_flank_chh_1500", 
"flagleaf_flank_chh_1600", "flagleaf_flank_chh_1700", "flagleaf_flank_chh_1800", 
"flagleaf_flank_chh_1900", "flagleaf_flank_chh_2000", "all3_avg_cg", 
"all3_avg_chg", "all3_avg_chh", "all3_flank_cg_100", "all3_flank_cg_200", 
"all3_flank_cg_300", "all3_flank_cg_400", "all3_flank_cg_500", 
"all3_flank_cg_600", "all3_flank_cg_700", "all3_flank_cg_800", 
"all3_flank_cg_900", "all3_flank_cg_1000", "all3_flank_cg_1100", 
"all3_flank_cg_1200", "all3_flank_cg_1300", "all3_flank_cg_1400", 
"all3_flank_cg_1500", "all3_flank_cg_1600", "all3_flank_cg_1700", 
"all3_flank_cg_1800", "all3_flank_cg_1900", "all3_flank_cg_2000", 
"all3_flank_chg_100", "all3_flank_chg_200", "all3_flank_chg_300", 
"all3_flank_chg_400", "all3_flank_chg_500", "all3_flank_chg_600", 
"all3_flank_chg_700", "all3_flank_chg_800", "all3_flank_chg_900", 
"all3_flank_chg_1000", "all3_flank_chg_1100", "all3_flank_chg_1200", 
"all3_flank_chg_1300", "all3_flank_chg_1400", "all3_flank_chg_1500", 
"all3_flank_chg_1600", "all3_flank_chg_1700", "all3_flank_chg_1800", 
"all3_flank_chg_1900", "all3_flank_chg_2000", "all3_flank_chh_100", 
"all3_flank_chh_200", "all3_flank_chh_300", "all3_flank_chh_400", 
"all3_flank_chh_500", "all3_flank_chh_600", "all3_flank_chh_700", 
"all3_flank_chh_800", "all3_flank_chh_900", "all3_flank_chh_1000", 
"all3_flank_chh_1100", "all3_flank_chh_1200", "all3_flank_chh_1300", 
"all3_flank_chh_1400", "all3_flank_chh_1500", "all3_flank_chh_1600", 
"all3_flank_chh_1700", "all3_flank_chh_1800", "all3_flank_chh_1900", 
"all3_flank_chh_2000", "SAM_avg_cg", "SAM_avg_chg", "SAM_avg_chh", 
"SAM_flank_cg_100", "SAM_flank_cg_200", "SAM_flank_cg_300", "SAM_flank_cg_400", 
"SAM_flank_cg_500", "SAM_flank_cg_600", "SAM_flank_cg_700", "SAM_flank_cg_800", 
"SAM_flank_cg_900", "SAM_flank_cg_1000", "SAM_flank_cg_1100", 
"SAM_flank_cg_1200", "SAM_flank_cg_1300", "SAM_flank_cg_1400", 
"SAM_flank_cg_1500", "SAM_flank_cg_1600", "SAM_flank_cg_1700", 
"SAM_flank_cg_1800", "SAM_flank_cg_1900", "SAM_flank_cg_2000", 
"SAM_flank_chg_100", "SAM_flank_chg_200", "SAM_flank_chg_300", 
"SAM_flank_chg_400", "SAM_flank_chg_500", "SAM_flank_chg_600", 
"SAM_flank_chg_700", "SAM_flank_chg_800", "SAM_flank_chg_900", 
"SAM_flank_chg_1000", "SAM_flank_chg_1100", "SAM_flank_chg_1200", 
"SAM_flank_chg_1300", "SAM_flank_chg_1400", "SAM_flank_chg_1500", 
"SAM_flank_chg_1600", "SAM_flank_chg_1700", "SAM_flank_chg_1800", 
"SAM_flank_chg_1900", "SAM_flank_chg_2000", "SAM_flank_chh_100", 
"SAM_flank_chh_200", "SAM_flank_chh_300", "SAM_flank_chh_400", 
"SAM_flank_chh_500", "SAM_flank_chh_600", "SAM_flank_chh_700", 
"SAM_flank_chh_800", "SAM_flank_chh_900", "SAM_flank_chh_1000", 
"SAM_flank_chh_1100", "SAM_flank_chh_1200", "SAM_flank_chh_1300", 
"SAM_flank_chh_1400", "SAM_flank_chh_1500", "SAM_flank_chh_1600", 
"SAM_flank_chh_1700", "SAM_flank_chh_1800", "SAM_flank_chh_1900", 
"SAM_flank_chh_2000", "gene_6_7_internode", "gene_7_8_internode", 
"gene_Vegetative_Meristem_16_19_Day", "gene_Ear_Primordium_2_4_mm", 
"gene_Ear_Primordium_6_8_mm", "gene_Embryo_20_DAP", "gene_Embryo_38_DAP", 
"gene_Endosperm_12_DAP", "gene_Endosperm_Crown_27_DAP", "gene_Germinatin_Kernels_2_DAI", 
"gene_Pericarp_Aleurone_27_DAP", "gene_Leaf_Zone_1__Symmetrical_", 
"gene_Leaf_Zone_2__Stomatal_", "gene_Leaf_Zone_3__Growth_", "gene_Mature_Leaf_8", 
"gene_Primary_Root_5_Days", "gene_Root___Cortex_5_Days", "gene_Root___Elongation_Zone_5_Days", 
"gene_Root___Meristem_Zone_5_Days", "gene_Secondary_Root_7_8_Days", 
"gene_B73_Mature_Pollen", "gene_Female_Spikelet_Collected_on_day_as_silk", 
"gene_Silk", "gene_coefvar", "gene_median", "gene_tau", "syntenicgene_6_7_internode", 
"syntenicgene_7_8_internode", "syntenicgene_Vegetative_Meristem_16_19_Day", 
"syntenicgene_Ear_Primordium_2_4_mm", "syntenicgene_Ear_Primordium_6_8_mm", 
"syntenicgene_Embryo_20_DAP", "syntenicgene_Embryo_38_DAP", "syntenicgene_Endosperm_12_DAP", 
"syntenicgene_Endosperm_Crown_27_DAP", "syntenicgene_Germinatin_Kernels_2_DAI", 
"syntenicgene_Pericarp_Aleurone_27_DAP", "syntenicgene_Leaf_Zone_1__Symmetrical_", 
"syntenicgene_Leaf_Zone_2__Stomatal_", "syntenicgene_Leaf_Zone_3__Growth_", 
"syntenicgene_Mature_Leaf_8", "syntenicgene_Primary_Root_5_Days", 
"syntenicgene_Root___Cortex_5_Days", "syntenicgene_Root___Elongation_Zone_5_Days", 
"syntenicgene_Root___Meristem_Zone_5_Days", "syntenicgene_Secondary_Root_7_8_Days", 
"syntenicgene_B73_Mature_Pollen", "syntenicgene_Female_Spikelet_Collected_on_day_as_silk", 
"syntenicgene_Silk", "syntenicgene_coefvar", "syntenicgene_median", 
"syntenicgene_tau", "TEfam_ear_primordium.2mm", "TEfam_ear_primordium.6mm", 
"TEfam_embryo_20d", "TEfam_embryo_38d", "TEfam_endosperm_12d", 
"TEfam_endosperm_crown", "TEfam_germinating.kernels_2d", "TEfam_internode_6to7", 
"TEfam_internode_7to8", "TEfam_leaf_8", "TEfam_leaf_growth.zone", 
"TEfam_leaf_stomatal.zone", "TEfam_leaf_symmetrical.zone", "TEfam_meristem_vegetative", 
"TEfam_pollen_mature", "TEfam_root_cortex", "TEfam_root_elongation.zone", 
"TEfam_root_meristem.zone", "TEfam_root_primary", "TEfam_root_secondary", 
"TEfam_seed_pericarp.aleurone", "TEfam_silk_mature", "TEfam_spikelet_female", 
"TEfamMedian", "TEfamMedianPerCopy", "TEfamMedianPerBp", "TEfam_tau", 
"cmmb", "subgenome", "famsize")
                              
categories$nicefeature=mapvalues(categories$feature, from=names(nicefeaturenames), to=nicefeaturenames)
         
library(stargazer)
stargazer(categories, summary=F, rownames=F, align=T)                                 

imp$category=mapvalues(imp$feat, from=categories$feature, to=categories$nicecategory)
imp$feat=mapvalues(imp$feat, from=as.character(categories$feature), to=as.character(categories$nicefeature))

#meltimp=melt(imp[rev(order(abs(imp$X.IncMSE))),][1:30,]) ## keep this managable!
#meltimpsum=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="abs") %>% summarize_if(.predicate="is.numeric", .funs="sum") %>% arrange(desc(X.IncMSE))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
#meltimpall=melt(imp[rev(order(abs(imp$X.IncMSE))),]) ## keep this managable!
meltimp=melt(imp[rev(order(abs(imp$scaled))),-'X.IncMSE'][1:30,]) ## keep this managable!
meltimpsum=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="sum", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
meltimpall=melt(imp[rev(order(abs(imp$scaled))),]) ## keep this managable!
meltimpmean=melt(imp %>% group_by(category) %>% summarize_if(.predicate="is.numeric", .funs="mean", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))


## for text:
library(stringr)
melt(imp %>% group_by(str_split_fixed(category, '_', 2)[,1]) %>% summarize_if(.predicate="is.numeric", .funs="mean", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))
melt(imp %>% group_by(str_split_fixed(category, '_', 2)[,1]) %>% summarize_if(.predicate="is.numeric", .funs="sum", na.rm=TRUE) %>% arrange(desc(scaled))) #sum=sum(X.IncMSE), meanscaled=mean(scaled)) %>% arrange(desc(sum))


## calc correlation with age for the 10 largest families in each sup??
feat30=data.frame(feat=as.character(imp$feat[rev(order(abs(imp$scaled)))[1:30]]))
### corr matrix by fam
for(fam in names(largest10)){
  feat30[,fam]<-NA
  for(feat in feat30$feat[-c(2,6)]){ ## have to get rid of sup and fam because they're not numeric
  feat30[feat30$feat==feat,fam]<-cor(ind$mya[ind$fam==fam], data.frame(ind)[,feat][ind$fam==fam]*1, use='na.or.complete')
    }
#  a=melt(round(cor(data.frame(ind)[ind$fam==fam, c('mya', feat30$feat)]*1, use='na.or.complete'),2))
#  
  }


## real plot
pdf(paste0('figure7.modeloutput.', Sys.Date(), '.pdf'), 30,12)

ispsc=ggplot(meltimp[meltimp$variable=='rmseMya',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimp$feat[order(meltimp[meltimp$variable=='rmseMya','value'])]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)

## make a correlation plot for alongside this one!!!! use feat30 i make above
mf=melt(feat30)
mf$sup=substr(mf$variable,1,3)
ispscCOR=ggplot(mf[mf$variable %in% names(largest5),], aes(x=factor(variable, levels=names(largest5)), y=factor(feat, levels=meltimp$feat[order(meltimp[meltimp$variable=='rmseMya','value'])]), size=abs(value), fill=factor(sign(value)))) +
                               geom_point(shape=21, stroke=0.1) + #scale_color_manual(values='black') + 
                              scale_fill_manual(name='Sign of Cor.', values=c('red', 'blue', 'black')) + 
                              ylab('') + xlab('') + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(), 
#                                                          axis.text.x=element_text(hjust=1, angle=90, size=rel(0.8))) +
                                                          axis.text.x=element_blank(), legend.position='none') +
                              labs(size=expression('r^2'))+
                              theme(axis.text.y=element_blank()) + 
                              annotate('rect', xmin = 0, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = dd.col['DHH'], alpha=0.3) +
                              annotate('rect', xmin = 5.5, xmax = 10.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTA'], alpha=0.3) +
                              annotate('rect', xmin = 10.5, xmax = 15.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTC'], alpha=0.3) +
                              annotate('rect', xmin = 15.5, xmax = 20.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTH'], alpha=0.3) +
                              annotate('rect', xmin = 20.5, xmax = 25.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTM'], alpha=0.3) +
                              annotate('rect', xmin = 25.5, xmax = 30.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTT'], alpha=0.3) +
                              annotate('rect', xmin = 30.5, xmax = 35.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTX'], alpha=0.3) +
                              annotate('rect', xmin = 35.5, xmax = 40.5, ymin = -Inf, ymax = Inf, fill = dd.col['RLC'], alpha=0.3) +
                              annotate('rect', xmin = 40.5, xmax = 45.5, ymin = -Inf, ymax = Inf, fill = dd.col['RLG'], alpha=0.3) +
                              annotate('rect', xmin = 45.5, xmax = 50.5, ymin = -Inf, ymax = Inf, fill = dd.col['RLX'], alpha=0.3) +
                              annotate('rect', xmin = 50.5, xmax = 55.5, ymin = -Inf, ymax = Inf, fill = dd.col['RIL'], alpha=0.3) +
                              annotate('rect', xmin = 55.5, xmax = 57.5, ymin = -Inf, ymax = Inf, fill = dd.col['RIT'], alpha=0.3) +
                              annotate('rect', xmin = 57.5, xmax = 62.5, ymin = -Inf, ymax = Inf, fill = dd.col['RST'], alpha=0.3)

ispscCOR4=ggplot(mf[mf$variable %in% names(largest5) & mf$sup %in% c('DHH', 'DTA', 'RLC', 'RLG'),], aes(x=factor(variable, levels=names(largest5)), y=factor(feat, levels=meltimp$feat[order(meltimp[meltimp$variable=='rmseMya','value'])]), size=abs(value), fill=factor(sign(value)))) +
                               geom_point(shape=21, stroke=0.01) + #scale_color_manual(values='black') + 
                              scale_fill_manual(name='Sign of Cor.', values=c('red', 'blue', 'black')) + 
                              ylab('') + xlab('') + theme(axis.title.x=element_blank(),axis.ticks.x=element_blank(), 
#                                                          axis.text.x=element_text(hjust=1, angle=90, size=rel(0.8))) +
                                                          axis.text.x=element_blank(), legend.position='none') +
                              labs(size=expression('r^2'))+
                              theme(axis.text.y=element_blank()) + 
                              annotate('rect', xmin = 0, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = dd.col['DHH'], alpha=0.3) +
                              annotate('rect', xmin = 5.5, xmax = 10.5, ymin = -Inf, ymax = Inf, fill = dd.col['DTA'], alpha=0.3) +
                              annotate('rect', xmin = 10.5, xmax = 15.5, ymin = -Inf, ymax = Inf, fill = dd.col['RLC'], alpha=0.3) +
                              annotate('rect', xmin = 15.5, xmax = 20.5, ymin = -Inf, ymax = Inf, fill = dd.col['RLG'], alpha=0.3)


#ispscCOR10=ggplot(mf[mf$variable %in% names(largest5),], aes(x=factor(variable, levels=names(largest5)), y=factor(feat, levels=rev(feat30$feat)), size=log10(abs(value)), color=sup, fill=factor(sign(value)))) +
#        geom_point(shape=21, stroke=2) + scale_color_manual(values=dd.col) + scale_fill_manual(values=c('red', 'blue', 'black'))
#
                      
                              
musc=ggplot(meltimpmean[meltimpmean$variable=='rmseMya',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpmean$category[order(abs(meltimpmean[meltimpmean$variable=='rmseMya','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)

## try out some model fits!
ss=fread('../age_model/segsites.bp.2019-01-31.txt')
mSSf=ggplot(ss, aes(segsites.bp, yhat.centered/2/3.3e-8/1e6, color=sup)) + 
                geom_line(aes(group = yhat.id), alpha = 0.1, data=ss[ss$yhat.id %in% sample(1:34151, 1000),]) + 
#                stat_summary(fun.y = mean, geom = "line",  size = 0.5, aes(group=fam, color=sup), alpha=0.2)+  
                stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  
                scale_color_manual(values=dd.col)+ theme(legend.position='none')
anth=fread('../age_model/anther_avg_chh.2019-02-01.txt')
mAnthf=ggplot(anth, aes(anther_avg_chh, yhat.centered/2/3.3e-8/1e6, color=sup)) + #
                geom_line(aes(group = yhat.id), alpha = 0.1, data=anth[anth$yhat.id %in% sample(1:34151, 1000),]) + 
#                stat_summary(fun.y = mean, geom = "line",  size = 0.5, aes(group=fam, color=sup), alpha=0.2)+  
                stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  
                scale_color_manual(values=dd.col) + theme(legend.position='none')
bp=fread('../figures/tebp.2019-07-15.txt') ## ugh put in the wrong directory :(
mTEBPf=ggplot(bp, aes(tebp, yhat.centered/2/3.3e-8/1e6, color=sup)) + #
                geom_line(aes(group = yhat.id), alpha = 0.1, data=bp[bp$yhat.id %in% sample(1:34151, 1000),]) + 
#                stat_summary(fun.y = mean, geom = "line",  size = 0.5, aes(group=fam, color=sup), alpha=0.2)+  
                stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  
                scale_color_manual(values=dd.col) + theme(legend.position='none')
close=fread('../figures/closest.2019-07-15.txt') ## ugh put in the wrong directory :(
mClosestf=ggplot(close, aes(closest, yhat.centered/2/3.3e-8/1e6, color=sup)) + #
                geom_line(aes(group = yhat.id), alpha = 0.1, data=close[close$yhat.id %in% sample(1:34151, 1000),]) + 
#                stat_summary(fun.y = mean, geom = "line",  size = 0.5, aes(group=fam, color=sup), alpha=0.2)+  
                stat_summary(fun.y = mean, geom = "line",  size = 2, aes(group=sup, color=sup))+  
                scale_color_manual(values=dd.col) + theme(legend.position='none')

                              
                              
                              
## plot the raw correlations
rSS=ggplot(ind, aes(x=segsites.bp, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth(se=F) + ylim(0,1) + theme(legend.position='none') + scale_color_manual(values=dd.col) + xlim(0,0.05)
rD=ggplot(ind, aes(group=disruptor, y=mya, color=sup)) + geom_boxplot() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)
rAnth=ggplot(ind, aes(x=anther_avg_chh, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth(se=F) + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)+ xlim(0,0.05)
rP=ggplot(ind, aes(x=TEfam_pollen_mature, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth(se=F) + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col) + xlim(0,2000)

rTEBP=ggplot(ind, aes(x=tebp, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth(se=F) + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)+ xlim(0,14308)
rClosest=ggplot(ind, aes(x=closest, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth(se=F) + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col) + xlim(0,112572)
                              

r=plot_grid(rSS + ylab('Age (Million years)') + xlab('Segregating sites per base pair'), 
            mSSf + ylim(0,5)+ ylab('Age (Million years)') + xlab('Permuted segregating sites per base pair'), 
            rAnth+ ylab('Age (Million years)') + xlab('CHH methylation (Anther tissue)'), 
            mAnthf + ylim(0,1)+ ylab('Age (Million years)') + xlab('Permuted CHH methylation (Anther tissue)'), 
            labels=c('D', 'E', 'F', 'G'), ncol=2, align='v')

r4=plot_grid(rSS + xlim(0,0.1)+ ylab('Age (Million years)') + xlab('Segregating sites per base pair'), 
             mSSf + ylim(0,5)+ xlim(0,0.1)+ ylab('Age (Million years)') + xlab('Permuted segregating sites per base pair'), 
             rAnth + xlim(0,0.05)+ ylab('Age (Million years)') + xlab('CHH methylation (Anther tissue)'), 
             mAnthf + ylim(0,1) + xlim(0,0.05)+ ylab('Age (Million years)') + xlab('Permuted CHH methylation (Anther tissue)'), 
             labels=c('D', 'E', 'F', 'G'), ncol=2, align='v')

rSupp=plot_grid(rTEBP+ ylab('Age (Million years)') + xlab('TE length (base pairs)'), 
                mTEBPf+ ylab('Age (Million years)') + xlab('Permuted TE length (base pairs)'), 
                rClosest+ ylab('Age (Million years)') + xlab('Distance to gene (base pairs)'),
                mClosestf+ ylab('Age (Million years)') + xlab('Permuted distance to gene (base pairs)'), 
                labels='AUTO', ncol=2, align='v')
### overleaf has a file size limit!
#plot_grid(musc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
#          ispsc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
#          r, 
#          ncol = 3, labels=c('A', 'B', ''), align = 'v', rel_widths=c(1,1,1.5))
    
plot_grid(musc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
          ispsc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
          ispscCOR,
          r4, 
          ncol = 4, labels=c('A', 'B', 'C', ''), align = 'v', rel_widths=c(1,1,2,1.75))
                                  
                              
## supp
#plot_grid(ispbysupscxM+ theme(legend.position='none'), ispbysupscxN+ theme(legend.position='none'), labels = "AUTO", ncol = 2, align = 'v', rel_widths=c(0.3,1))

dev.off()

## make a reasonably sized png
png(paste0('figure6.modeloutput.', Sys.Date(), '.png'), 30, 12, units='in', res=600)#*300,12*300) ## *300 dpi

                              
plot_grid(musc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
          ispsc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
          ispscCOR,
          r4, 
          ncol = 4, labels=c('A', 'B', 'C', ''), align = 'h', axis='t', rel_widths=c(1,1,2,1.75), scale=c(1,1,1,1))
dev.off()
                              
## make a reasonably sized png
png(paste0('figure6.modeloutput4sups.', Sys.Date(), '.png'), 30, 12, units='in', res=600)#*300,12*300) ## *300 dpi

plot_grid(musc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
          ispsc + theme(legend.position='none') + ylab('Reduction in square root mean squared error (Mya)') + xlab(''), 
          ispscCOR4,
          r4, 
          ncol = 4, labels=c('A', 'B', 'C', ''), align = 'h', axis='t',rel_widths=c(1,1,2,1.75), scale=c(1,1,1,1))
dev.off()
  
## make a reasonably sized png
png(paste0('supp.modeloutput.', Sys.Date(), '.png'), 30, 12, units='in', res=600)#*300,12*300) ## *300 dpi

rSupp
dev.off()

## old stuff!
############ start plotting
##pdf(paste0('figure6.modeloutput.', Sys.Date(), '.pdf'), 30,12)
ip=ggplot(impsort, aes(y=category, x=sum)) + geom_point() + scale_y_discrete(limits=impsort$category)
#isp=ggplot(impsort[impsort$meanscaled>1.9,], aes(y=category, x=meanscaled)) + geom_point() + 
#       scale_y_discrete(limits=rev(impsort[impsort$meanscaled>1.9,] %>% arrange(desc(meanscaled)))$category)
isp=ggplot(meltimp[meltimp$variable=='X.IncMSE',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimp$feat[order(abs(meltimp[meltimp$variable=='X.IncMSE','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)
su=ggplot(meltimpsum[meltimpsum$variable=='X.IncMSE',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpsum$category[order(abs(meltimpsum[meltimpsum$variable=='X.IncMSE','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)
#su=ggplot(impsort[impsort$meanscaled>1.9,], aes(y=

ispbysup=ggplot(meltimpall[meltimpall$variable=='X.IncMSE',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimp$feat[order(meltimp[meltimp$variable=='X.IncMSE','value'])]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free')
       
       
       
       
ispsc=ggplot(meltimp[meltimp$variable=='scaled',], aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimp$feat[order(meltimp[meltimp$variable=='scaled','value'])]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)
susc=ggplot(meltimpsum[meltimpsum$variable=='scaled',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpsum$category[order(abs(meltimpsum[meltimpsum$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)

musc=ggplot(meltimpmean[meltimpmean$variable=='scaled',], aes(x=category, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
       scale_x_discrete(limits=meltimpmean$category[order(abs(meltimpmean[meltimpmean$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors)


ispbysupsc=ggplot(meltimpall[meltimpall$variable=='scaled',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free')


ispbysupscx=ggplot(meltimpall[meltimpall$variable=='scaled',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free_y')
ispbysupscxM=ggplot(meltimpall[meltimpall$variable=='scaled' & meltimpall$category=='flank_methylation_mnase',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free_y')
ispbysupscxN=ggplot(meltimpall[meltimpall$variable=='scaled' & meltimpall$category!='flank_methylation_mnase',][order(abs(meltimpall[meltimpall$variable=='scaled','value'])),], 
      aes(x=feat, y=abs(value), fill=category)) + geom_bar(position="dodge",stat="identity") + 
#       geom_errorbar(aes(ymin=weight-std/2, ymax=weight+std/2), size=.3, width=.2, position=position_dodge(.9)) +
#       scale_x_discrete(limits=meltimpall$feat[order(abs(meltimpall[meltimpall$variable=='scaled','value']))]) +  
       coord_flip() + #scale_fill_brewer(palette='Set3') #+ colScale
       scale_fill_manual(values=myColors) + facet_wrap(~category, scales='free_y')


plot_grid(su + theme(legend.position='none'), isp + theme(legend.position='none'), ispbysup + theme(legend.position='none'), labels = "AUTO", ncol = 3, align = 'v', rel_widths=c(1,1,2))

plot_grid(susc + theme(legend.position='none'), ispsc + theme(legend.position='none'), ispbysupsc + theme(legend.position='none'), labels = "AUTO", ncol = 3, align = 'v', rel_widths=c(1,1,2))



## plot the raw correlations
rSS=ggplot(ind, aes(x=segsites.bp, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none') + scale_color_manual(values=dd.col) + xlim(0,0.05)
rD=ggplot(ind, aes(group=disruptor, y=mya, color=sup)) + geom_boxplot() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)
rAnth=ggplot(ind, aes(x=anther_avg_chh, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col)+ xlim(0,0.05)
rP=ggplot(ind, aes(x=TEfam_pollen_mature, y=mya, color=sup)) + geom_point(size=0.1, alpha=0.1) + stat_smooth() + ylim(0,1) + theme(legend.position='none')+ scale_color_manual(values=dd.col) + xlim(0,2000)

r=plot_grid(rSS, rD, rAnth, rP, ncol=2, align='v')

plot_grid(musc + theme(legend.position='none'), ispsc + theme(legend.position='none'), r, ncol = 3, align = 'v', rel_widths=c(1,1,1.5))
          
## supp
plot_grid(ispbysupscxM+ theme(legend.position='none'), ispbysupscxN+ theme(legend.position='none'), labels = "AUTO", ncol = 2, align = 'v', rel_widths=c(0.3,1))

dev.off()




### pollen supplement
po=ind %>% group_by(sup, fam,famsize) %>% dplyr::summarize(medianage=median(mya), pollen=median(TEfam_pollen_mature))

pdf('pollen_vs_age.pdf',14,8)
ggplot(po[po$famsize>=10 & !is.na(po$pollen),], aes(x=pollen, y=medianage, col=sup)) + geom_point()  + theme(legend.position='none') + scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F)
poall=ggplot(po[po$famsize>=10 & !is.na(po$pollen),], aes(x=pollen/famsize, y=medianage, col=sup)) + geom_point()  + theme(legend.position='none') + scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F)
poall
ggplot(po[po$famsize>=10 & !is.na(po$pollen),], aes(x=pollen, y=medianage, col=sup)) + geom_point()  + theme(legend.position='none') + scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F) + xlim(0,2) + ylim(0,2)
ggplot(po[po$famsize>=10 & !is.na(po$pollen),], aes(x=pollen/famsize, y=medianage, col=sup)) + geom_point()  + theme(legend.position='none') + scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F) + xlim(0,2) + ylim(0,2)
posmall=ggplot(po[po$famsize>=10 & !is.na(po$pollen),], aes(x=pollen/famsize, y=medianage, col=sup)) + geom_point()  + theme(legend.position='none') + scale_color_manual(values=dd.col) + geom_smooth(method='gam', se=F) + xlim(0,0.5) + ylim(0,2)
posmall
plot_grid(poall, posmall, labels='AUTO')                              
dev.off()


