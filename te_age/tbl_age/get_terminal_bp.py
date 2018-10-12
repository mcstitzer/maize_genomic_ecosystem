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

slen=int(sys.argv[2])

fams=[]
for te in b.keys():
	if te[0:8] not in fams:
		fams.append(te[0:8])
		if slen > 0:
			print '>'+te+'\n'+printFasta(''.join(b[te])[:slen])
		elif slen <=0:
			print '>'+te+'\n'+printFasta(''.join(b[te])[slen:])


#for entry in b:
#	if entry[0:8] in te:
#		print '>'+entry+'\n'+printFasta(''.join(b[entry]))
