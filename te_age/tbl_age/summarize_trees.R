library('ape')
library('dplyr')

rst=read.tree('RST.tre')
rit=read.tree('RIT.tre')
ril=read.tree('RIL.tre')
dta=read.tree('DTA.tre')
dtc=read.tree('DTC.tre')
dth=read.tree('DTH.tre')
dtm=read.tree('DTM.tre')
dtt=read.tree('DTT.tre')
dtx=read.tree('DTX.tre')
rlc=read.tree('RLC.1000.tre')
rlg=read.tree('RLG.1000.tre')
rlx=read.tree('RLX.1000.tre')
dhh=read.tree('DHH.last1000.tre')

## get the terminal branch length for each superfamily, put in one data frame
tes=data.frame(TEID=rst$tip.label, tbl=rst$edge.length[rst$edge[,2]<=Ntip(rst)])
tes=rbind(tes,data.frame(TEID=dta$tip.label, tbl=dta$edge.length[dta$edge[,2]<=Ntip(dta)]))
tes=rbind(tes,data.frame(TEID=dtc$tip.label, tbl=dtc$edge.length[dtc$edge[,2]<=Ntip(dtc)]))
tes=rbind(tes,data.frame(TEID=dth$tip.label, tbl=dth$edge.length[dth$edge[,2]<=Ntip(dth)]))
tes=rbind(tes,data.frame(TEID=dtm$tip.label, tbl=dtm$edge.length[dtm$edge[,2]<=Ntip(dtm)]))
tes=rbind(tes,data.frame(TEID=dtt$tip.label, tbl=dtt$edge.length[dtt$edge[,2]<=Ntip(dtt)]))
tes=rbind(tes,data.frame(TEID=dtx$tip.label, tbl=dtx$edge.length[dtx$edge[,2]<=Ntip(dtx)]))
tes=rbind(tes,data.frame(TEID=dhh$tip.label, tbl=dhh$edge.length[dhh$edge[,2]<=Ntip(dhh)]))                  
tes=rbind(tes,data.frame(TEID=rit$tip.label, tbl=rit$edge.length[rit$edge[,2]<=Ntip(rit)]))
tes=rbind(tes,data.frame(TEID=ril$tip.label, tbl=ril$edge.length[ril$edge[,2]<=Ntip(ril)]))
tes=rbind(tes,data.frame(TEID=rlc$tip.label, tbl=rlc$edge.length[rlc$edge[,2]<=Ntip(rlc)]))
tes=rbind(tes,data.frame(TEID=rlg$tip.label, tbl=rlg$edge.length[rlg$edge[,2]<=Ntip(rlg)]))          
tes=rbind(tes,data.frame(TEID=rlx$tip.label, tbl=rlx$edge.length[rlx$edge[,2]<=Ntip(rlx)]))          

     
## get rid of reversed notation from mafft
tes$TEID=gsub('_R_','',tes$TEID)
tes=tes[order(tes$TEID),]
write.table(tes, 'B73v4_terminalbranchlength.txt', col.names=T, row.names=F, sep='\t', quote=F)


tes$sup=substr(tes$TEID,1,3)
tes$fam=substr(tes$TEID,1,8)
te_fam=tes %>% group_by(sup,fam) %>% dplyr::summarize(avg_tbl=mean(tbl, na.rm=T))
write.table(te_fam, 'B73v4_terminalbranchlength.fam.txt', quote=F, sep='\t', col.names=T, row.names=F)

## note i'm taking averages of averages here! Should each family count equally? Or should it be weighted by famsize?
te_sup=te_fam %>% group_by(sup) %>% dplyr::summarize(avg_tbl=mean(avg_tbl, na.rm=T))
write.table(te_sup, 'B73v4_terminalbranchlength.sup.txt', quote=F, sep='\t', row.names=F, col.names=T)


tes$famsize=table(tes$fam)[tes$fam]
