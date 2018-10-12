import sys

#python split_fasta.py TAIR10_Chr.all.allsubtracted.ltrdigest_5ltr.newName.fa TAIR10_Chr.all.allsubtracted.ltrdigest_3ltr.newName.fa

fasta=open(sys.argv[1], 'r')
out=None
for line in fasta:
	if line.startswith('>'):
		if out==None:
			out=open(sys.argv[2]+'/'+line.strip()[1:]+'.fa', 'w')
		else:
			out.close()
			out=open(sys.argv[2]+'/'+line.strip()[1:]+'.fa', 'w')
		out.write(line)
	else:
		out.write(line)
out.close()
