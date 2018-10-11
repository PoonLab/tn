#Usage "python3 fasTOcsv.py3 input.fas output.csv"

import sys, os
infile = open(sys.argv[1], "r")
outfile = open('csvFromAlignedFas.csv', "a+")

outtext = 'id,year,sequence'

for line in infile:
	if line.startswith(">"):	
		seq_id = line.strip('>').split('_')[0]
		seq_date = line.strip().split('_')[1]		
		outtext += ('\n' + seq_id + ',' + seq_date + ',')
	else:
		outtext += (line.strip())

outfile.write(outtext)
outfile.close()

