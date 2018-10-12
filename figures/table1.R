library(stargazer)
library(data.table)
library(dplyr)

source('../GenomeInfo.R') ## we'll need TESUPFACTORLEVELS here!
source('color_palette.R')

## note that this expects there to be an ind data frame

## could make these use GENOME to get the file name
techar=fread('../te_characteristics/B73_TE_individual_copies.2018-09-19.txt')


d = techar %>% group_by(sup) %>% dplyr::summarize(number.copies=n(), 
                                      number.families=length(unique(fam))
                                      )


d=data.frame(class=c('DNA transposon', 'DNA transposon','DNA transposon','DNA transposon','DNA transposon','DNA transposon','DNA transposon','Retrotransposon','Retrotransposon','Retrotransposon','Retrotransposon','Retrotransposon','Retrotransposon'), 
             order=c('Helitron', 'TIR', 'TIR', 'TIR', 'TIR', 'TIR', 'TIR', 'nonLTR', 'nonLTR', 'LTR', 'LTR', 'LTR', 'nonLTR'),
             common.name=c('Helitron', 'hAT', 'CACTA', 'Pif/Harbinger', 'Mutator', 'Tc1/Mariner', 'Unknown TIR', 'Jockey', 'L1', 'Copia', 'Gypsy', 'Unknown LTR', 'SINE'), 
             d)
d=d[match(TESUPFACTORLEVELS, d$sup),]             
d=d[,c('class', 'order', 'sup', 'common.name', 'number.copies', 'number.families')]
names(d)=c('Class', 'Order', 'Superfamily', 'Common Name', 'Number Copies', 'Number Families')
stargazer(d, summary=F, rownames=F, align=T)


## then in doc, I add caption{Superfamilies in the maize genome} and after label, use \scalebox{0.5}{, which ends after \end{tabular} with a } before \end{table}
