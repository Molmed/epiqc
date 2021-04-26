#! /usr/bin/env python3
# --- combine methylation output files ------------------------------------------------------------ #
#
#  e.g., for file 1:
#  chr1	10649	10670	75.00	3	1
#  and file 2:
#  chr1	10649	10670	50.00	2	2
#  will produce:
#  chr1	10649	10670	62.50	5	3
#  note that output is random, might need to sort
#
#  usage:  python combine_methylation_outputs.py <sample_outname> <file1> <file2> ... <fileN>
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
			CpG_values = (int(line[4]), int(line[5]))
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
		with open(biscov, 'r') as infile:
			merger(infile)

# output
with gzip.open(sample_outname+'.bedGraph.gz', 'wt') as outfile:
	for locus in combined_data:
		try:
			pct = round(float(combined_data[locus][0]) / (float(combined_data[locus][0]) + float(combined_data[locus][1])) * 100.0, 3)
			outfile.write(locus.replace('.','\t') +'\t'+ str(pct) +'\t'+ str(combined_data[locus][0]) +'\t'+ str(combined_data[locus][1]) +'\n')
		except ZeroDivisionError:
			pass
		

