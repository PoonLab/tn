#Working towards a simple method of fasta header pre-processing (purely for

#USAGE: "python3 prePro.py3 [input file]"

import sys, os
inFile = open(sys.argv[1], "r")
#idLoc = sys.argv[2]
#cyLoc = sys.argv[3]
outfile = open(sys.argv[2]+"PP.fas", "a+")

#print(idLoc)
#print(cyLoc)

for line in inFile: 
	if line.startswith(">"):
		line = line.strip('>').strip('\n')

		seqId = line.split('_')[1]
		#for i in range(0,int(len(idLoc)/2)):
			#seqId = seqId.split(idLoc[i*2])[int(idLoc[i*2+1])]

		seqCy = line.split('_')[len(line.split('_'))-1]
		#for i in range(0,int(len(cyLoc)/2)):
			#seqCy = seqCy.split(cyLoc[i*2])[int(cyLoc[i*2+1])]


		line = (">" + seqId + "_" + str(seqCy) + '\n')
		outfile.write(line)
		
	else: 
		outfile.write(line)
