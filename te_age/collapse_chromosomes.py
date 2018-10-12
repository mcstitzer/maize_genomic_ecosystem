import sys
import textwrap

def readFasta(filename):
	fastafile=open(filename, 'r')
	fastadict={}
	for line in fastafile:		
		if line.startswith('>') and line.strip()[1:] not in fastadict:
			seqname=line.strip()[1:]
			fastadict[seqname]=[]
		elif line.startswith('>') and line.strip()[1:] in fastadict:
			seqname=line.strip()[1:]
#			continue
		else:
			fastadict[seqname].append(line.strip())
	for entry in fastadict:
		fastadict[entry]=''.join(fastadict[entry])
	return fastadict

def printFasta(sequence, width=70):
	return '\n'.join(sequence[i:i+width] for i in range(0,len(sequence), width))


b=readFasta(sys.argv[1])
#out=open('stupid_test.txt', 'w')
#out.write('read in all the fasta entries')
for entry in b:
	print '>'+entry+'\n'+printFasta(''.join(b[entry]))


superfams=['DHH', 'DTA',  'DTC',  'DTH',  'DTM',  'DTT',  'DTX',  'RIL', 'RIT',  'RLC',  'RLG',  'RLX',  'RST']

DHH=open('tbl_age/DHH.fa', 'w')
DTA=open('tbl_age/DTA.fa', 'w')
DTC=open('tbl_age/DTC.fa', 'w')
DTH=open('tbl_age/DTH.fa', 'w')
DTM=open('tbl_age/DTM.fa', 'w')
DTT=open('tbl_age/DTT.fa', 'w')
DTX=open('tbl_age/DTX.fa', 'w')
RIL=open('tbl_age/RIL.fa', 'w')
RIT=open('tbl_age/RIT.fa', 'w')
RLC=open('tbl_age/RLC.fa', 'w')
RLG=open('tbl_age/RLG.fa', 'w')
RLX=open('tbl_age/RLX.fa', 'w')
RST=open('tbl_age/RST.fa', 'w')


for entry in b:
        if entry[0:3]=='DHH':
                DHH.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='DTA':
                DTA.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='DTC':
                DTC.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='DTH':
                DTH.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='DTM':
                DTM.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='DTT':
                DTT.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='DTX':
                DTX.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='RIL':
                RIL.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='RIT':
                RIT.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='RLC':
                RLC.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='RLG':
                RLG.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='RLX':
                RLX.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')
        if entry[0:3]=='RST':
                RST.write('>'+entry+'\n'+printFasta(''.join(b[entry]))+'\n')


DHH.close()
DTA.close()
DTC.close()
DTH.close()
DTM.close()
DTT.close()
DTX.close()
RIL.close()
RIT.close()
RLC.close()
RLG.close()
RLX.close()
RST.close()
