#! /usr/bin/env python
#  METHylation TRANSlation:
#  convert BSSEEKER2 to BEDGRAPH

# IN:
# chr1    C       10497   CG      CG      1.0     1       1
# chr1    C       10502   CHG     CT      0.0     0       1
# chr1    C       10506   CHH     CC      0.0     0       1
# chr1    C       10507   CHG     CT      0.0     0       1
# chr1    C       10517   CHG     CT      0.0     0       2
# chr1    C       10522   CHH     CT      0.0     0       2
# chr1    C       10524   CHG     CC      0.0     0       2
# chr1    C       10525   CG      CG      1.0     2       2
# OUT:
# chr1	10497	10497   100.0   1   0
# chr1	10525	10525	100.0	2	0

import os
import sys
import gzip

fullPath = sys.argv[1]
bsseekerFile = os.path.basename(fullPath)
direc = os.path.dirname(fullPath)

with gzip.open(direc+'/'+bsseekerFile.replace('CGmap','bedGraph'), 'wt') as outfile:
	with gzip.open(fullPath, 'rt') as infile:
		for line in infile:
			line = line.strip().split('\t')
			if line[3] == "CG":
				chr = line[0]
				pos_end = line[2]
				pos_start = str(int(line[2])-1)
				pct = str(float(line[5])*100.0)
				freqC = line[6]
				freqT = int(line[7]) - int(freqC)
				outfile.write(chr + '\t' + pos_start + '\t' + pos_end +'\t' + pct + '\t' + str(freqC) + '\t' + str(freqT) + '\n')
