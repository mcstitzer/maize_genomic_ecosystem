## Mnase data from Rogers-Melnick et al 2016

This is a data set where I really should remap for biological interpretation, but I'm not finding the code to do so in a quick search. For now, I'll convert coordinates, but consider switching in the future.

1. get data: 
 - Shoot ```wget http://de.iplantcollaborative.org/dl/d/E9018AC6-8460-4374-AF1C-4EABF92F1295/AP.bfthresh1.1.MNaseHS.Ranges.dat```
 - Root ```wget http://de.iplantcollaborative.org/dl/d/4517FAAA-8969-4B43-AE2E-E1A91DD65BFB/RP.bfthresh1.1.MNaseHS.Ranges.dat```
 Note that these are chromosome, start, end, 0-based end exclusive
 Also, these are the MNase hotspots, so I'll just look for overlaps with interesting features

2. convert to AGPv4 coordinates on [ensembl](http://plants.ensembl.org/Zea_mays/Tools/AssemblyConverter)
 - ensembl doesn't like these files end in .dat, so append a .bed first. 
```mv AP.bfthresh1.1.MNaseHS.Ranges.dat AP.bfthresh1.1.MNaseHS.Ranges.dat.bed``` 
```mv RP.bfthresh1.1.MNaseHS.Ranges.dat RP.bfthresh1.1.MNaseHS.Ranges.dat.bed```

## the output in v4 coordinates are ```AP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed``` and ```RP.bfthresh1.1.MNaseHS.Ranges.AGPv4.bed```
I included these files in this repo because I am not sure of the reproducibility of the ensembl conversion

