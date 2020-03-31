import os.path, re, glob

array_dir = '../new_array_data'
bed_dir = '../new_bedgraph'


bedgraph_files = glob.glob(os.path.join(bed_dir, '*.bedGraph.gz'))

strandmerged = ['combined/strandmerged_' + os.path.basename(x).split('.')[0] + '.cov.gz' for x in bedgraph_files]

reps_attr = list(set(['_'.join(x[:3]) for x in [os.path.basename(x).split('_') for x in bedgraph_files]]))

reps_combined = ['combined/repsCombined_' + x + '.cov.gz' for x in reps_attr]

sample_attr = list(set(['_'.join([x[0], 'HG001-7', x[2]]) for x in [x.split('_') for x in reps_attr]]))

samples_combined = ['combined/samplesCombined_' + x + '.cov.gz' for x in sample_attr]

localrules: all

rule all:
	input:
		strandmerged,
		reps_combined,
		samples_combined,
		'array_anot/array_anot_with_hg38.Rd',
		'bsseq_robjects/epiqc_bsseq_array_sites_reps_merged.Rd',
		'bsseq_robjects/mapping_table_epic_array_hg38_hg19.txt'

rule cg_coords:
	output:
		'cg_coords/cg_coords.Rd'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {output} < cg_coords.R
		"""


rule sum_dinuc:
	input:
		os.path.join(bed_dir, '{bedgraph}.bedGraph.gz'),
		'cg_coords/cg_coords.Rd'
	output:
		'combined/strandmerged_{bedgraph}.cov.gz'

	shell:
		"""
		module load R
		R --no-restore --no-save --args {input[1]} {input[0]} {output} < strand_merge.R
		"""

def get_strandmerged(wildcards):
	files = strandmerged
	patt = re.compile(wildcards.prep + '_' + wildcards.sample + '_' + wildcards.lab)
	inp = []
	for f in files:
		mobj = patt.search(f)
		if mobj:
			inp.append(f)
	return(inp)


rule combine_reps:
	input:
		get_strandmerged

	output:
		'combined/repsCombined_{prep}_{sample}_{lab}.cov.gz'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < combine_reps.R
		"""

def get_combined(wildcards):
	files = reps_combined
	patt = re.compile(wildcards.prep + '_.+_' + wildcards.lab)
	inp = []
	for f in files:
		mobj = patt.search(f)
		if mobj:
			inp.append(f)
	return(inp)

rule combine_samples:
	input:
		get_combined

	output:
		'combined/samplesCombined_{prep}_HG001-7_{lab}.cov.gz'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < combine_reps.R
		"""

rule liftover_array:
	input:
		os.path.join(array_dir, 'rgset.RData'),
		os.path.join(array_dir, 'hg19ToHg38.over.chain')

	output:
		'array_anot/array_anot_with_hg38.Rd'

	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < liftover_200327.R
		"""

rule bsseq_se_array_sites:
	input:
		reps_combined,
		'array_anot/array_anot_with_hg38.Rd'

	output:
		'bsseq_robjects/epiqc_bsseq_array_sites_reps_merged.Rd',
		'bsseq_robjects/mapping_table_epic_array_hg38_hg19.txt'

	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < bsseq_array_sites.R
		"""	
