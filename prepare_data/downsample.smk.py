import os.path, glob, decimal

covgs1 = ['5', '10', '20', '30', '40', '50']
covgs2 = ['5', '10', '20', '30', '40', '50', '100', '150', '200']

bed_dir = 'combined/'

localrules: all, cat_coverage_table, cat_coverage_table_post, cat_coverage_table_cumulative

def split_file_names(flist):
	spl =  [os.path.basename(x).split('.')[0].split('_') for x in flist]
	attr = ['_'.join([x[1], x[2], x[3]]) for x in spl]
	return({k:v for (k, v) in zip(attr, flist)})

reps_combined = glob.glob(os.path.join(bed_dir, 'repsCombined_*.cov.gz'))
reps_attr = split_file_names(reps_combined)
samples_combined = glob.glob(os.path.join(bed_dir, 'samplesCombined_*.cov.gz'))
samples_attr = split_file_names(samples_combined)
attr = {**reps_attr, **samples_attr}

prep_sample_lab_cov = [x + '_' + y for x in reps_attr.keys() for y in covgs1] + [x + '_' + y for x in samples_attr.keys() for y in covgs2]

sample_lab_cov = list(set(['_'.join([y[1], y[2], y[3]]) for y in [x.split('_') for x in prep_sample_lab_cov]]))

rule all:
	input:
		['downsampled/downsampled_' + x + 'x.cov.gz' for x in prep_sample_lab_cov],
		'total_coverage/coverage_table.txt',
		'cumulative_coverage/cumulative_coverage_table.txt',
		'total_coverage/coverage_table_post.txt' ,
		['upset_matrices/upset_' + x + 'x.txt.gz' for x in sample_lab_cov]



rule determine_coverage:
	input:
		lambda wildcards: attr[wildcards.sample_key]
	output:
		temp('total_coverage/{sample_key}_total_coverage.txt')
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < file_coverage.R
		"""

rule cat_coverage_table:
	input:
		expand('total_coverage/{sample_key}_total_coverage.txt', sample_key=attr.keys())
	output:
		'total_coverage/coverage_table.txt',
	shell:
		"""
		cat {input} > {output[0]}
		"""

def get_coverage1(wildcards, cov_file):
	decimal.getcontext().prec = 8
	covs = []
	sample_coverage = None
	for line in open(cov_file, 'r'):
		prep = line.split()[0]
		sample = line.split()[1]
		lab = line.split()[2]
		coverage = line.split()[3]
		if sample == wildcards.sample and lab == wildcards.lab:
			covs.append(int(coverage))
			if prep in wildcards.prep:
				sample_coverage = coverage

	min_cov = min(covs)
	fraction = str(decimal.Decimal(min_cov)/decimal.Decimal(sample_coverage))
	return(fraction)

def get_coverage2(wildcards, cov_file):
	decimal.getcontext().prec = 8
	print(wildcards)
	print(cov_file)
	sample_coverage = None
	for line in open(cov_file, 'r'):
		prep = line.split()[0]
		sample = line.split()[1]
		lab = line.split()[2]
		coverage = line.split()[4]
		if sample == wildcards.sample and prep == wildcards.prep and lab == wildcards.lab:
				sample_coverage = coverage
				print([sample, prep, lab])
				break
	print(sample_coverage)
	fraction = str(decimal.Decimal(wildcards.cov)/decimal.Decimal(sample_coverage))
	print(fraction)
	return(fraction)

rule downsample_bedgraph:
	input:
		lambda wildcards: attr['_'.join([wildcards.prep, wildcards.sample, wildcards.lab])],
		'total_coverage/coverage_table.txt'

	params:
		sample_fraction = lambda wildcards, input: get_coverage2(wildcards, input[1])

	output:
		'downsampled/downsampled_{prep}_{sample}_{lab}_{cov}x.cov.gz'

	run:
		if params.sample_fraction < 1:
			shell('zcat {input[0]} | python3 downsample_methylKit.py --fraction {params.sample_fraction} --bedGraph | gzip > {output}')
		else:
			shell('touch {output}')
		


rule determine_coverage_post:
	input:
		'downsampled/downsampled_{prep}_{sample}_{lab}_{cov}x.cov.gz'
	output:
		temp('total_coverage/{prep}_{sample}_{lab}_{cov}_total_coverage.txt')
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < file_coverage.R
		"""

rule cat_coverage_table_post:
	input:
		['total_coverage/' + x + '_total_coverage.txt' for x in prep_sample_lab_cov],

	output:
		'total_coverage/coverage_table_post.txt',

	shell:
		"""
		cat {input} > {output[0]}
		"""

rule cumulative_coverage:
	input:
		'downsampled/downsampled_{prep}_{sample}_{lab}_{cov}x.cov.gz'
	output:
		temp('cumulative_coverage/{prep}_{sample}_{lab}_{cov}x_cumulative_coverage.txt')
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < coverage_table.R
		"""

rule cat_coverage_table_cumulative:
	input:
		['cumulative_coverage/' + x + 'x_cumulative_coverage.txt' for x in prep_sample_lab_cov],

	output:
		'cumulative_coverage/cumulative_coverage_table.txt',
	shell:
		"""
		cat {input} > {output[0]}
		"""

rule cg_coords:
	output:
		'cg_coords/cg_coords.Rd'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {output} < cg_coords.R
		"""

def upset_input(wildcards):
	inp = []
	patt = re.compile(wildcards.slc + '$')
	for pslc in prep_sample_lab_cov:
		mobj = patt.search(pslc)
		if mobj:
			inp.append('downsampled/downsampled_' + pslc + 'x.cov.gz')
	return(inp)


rule upset_matrix:
	input:
		upset_input,
		'cg_coords/cg_coords.Rd'
	output:
		'upset_matrices/upset_{slc}x.txt.gz'

	params:
		cov = lambda wildcards: wildcards.slc.split('_')[2]

	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {params.cov} {output} < upset_matrix.R
		"""
