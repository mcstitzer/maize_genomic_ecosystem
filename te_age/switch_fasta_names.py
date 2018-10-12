import sys


fa=open(sys.argv[1], 'r')
tab=open(sys.argv[2], 'r')

origname=[]
newname=[]

for line in tab:
	fields=line.strip().split()
	origname.append(fields[0])
	newname.append(fields[1])
tab.close()

use=0
for line in fa:
	if line.startswith('>'):
		name=line.strip()[1:]    ## get just the name
		if name in origname:
			print ">"+newname[origname.index(name)]
			use=1
		else:
			use=0
	elif use==1:
		print line.strip()
	elif use==0:
		continue
