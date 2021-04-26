#! /usr/bin/env python
#	1/27/20

# --- combine methylation output files ------------------------------------------------------------ #
#
#	e.g., for file 1:
#	chr13	48423360	48423361	.	1	-	48423360	48423361	0,0,0	1	100
#	and file 2:
#	chr13	48423360	48423361	.	2	-	48423360	48423361	0,0,0	2	50
#	will produce:
#	chr13	48423360	48423361	66.67	2	1
#	this is generally useful for ONT data processed by e.g. Megalodon which reports in methylBedGraph
#
#	usage:	python methylbedGraph_combine.py <sample_outname> <file1> <file2> ... <fileN>
#
# ------------------------------------------------------------------------------------------------- #

import gzip
import sys

sample_outname = sys.argv[1]
biscov_files   = sys.argv[2:]

combined_data = {}

def merger(biscov):
	for line in infile:
		if not 'track' in line:
			line = line.strip().split('\t')
			locus = line[0]+'.'+line[1]+'.'+line[2]
	
		# convert # bases and pct --> Cs and Ts
		numBases = float(line[9])
		if numBases > 0:
			pct = float(line[10]) / 100.0
			numCs = int(round(float(pct * numBases)))
			numTs = int(numBases - numCs)
			CpG_values = (numCs, numTs)
			if locus not in combined_data:
				combined_data[locus] = CpG_values
			else:
				combined_data[locus] = tuple(map(sum, zip(combined_data[locus], CpG_values)))

# process
for biscov in biscov_files:
	print('incorporating ' + biscov + ' ...')
	if biscov.endswith('.gz'):
		with gzip.open(biscov, 'rt') as infile:
			merger(infile)
	else:
		with open(biscov, 'rt') as infile:
			merger(infile)

# output
print('writing output ...')
with gzip.open(sample_outname+'.bedGraph.gz', 'wt') as outfile:
	for locus in combined_data:
		pct = round(float(combined_data[locus][0]) / (float(combined_data[locus][0]) + float(combined_data[locus][1])) * 100.0, 3)
		outfile.write(locus.replace('.','\t') +'\t'+ str(pct) +'\t'+ str(combined_data[locus][0]) +'\t'+ str(combined_data[locus][1]) +'\n')

