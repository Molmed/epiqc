import os.path, glob

samples1 = ['HG00' + str(n) for n in range(1,8)]
samples2 = 'HG001-7'
preps = ['EMSeq', 'MethylSeq', 'SPLAT', 'TruSeq']
covgs1 = ['5', '10', '15', '20']
covgs2 = ['25', '50', '80', '100', '140']

bed_dir = './combined/'

merged_reps = [os.path.basename(x).split('.')[0] for x in  glob.glob(os.path.join(bed_dir, 'repsCombined_*.cov.gz'))]
merged_samples = [os.path.basename(x).split('.')[0] for x in  glob.glob(os.path.join(bed_dir, 'samplesCombined_*.cov.gz'))]


localrules: all, cat_coverage_table, cat_coverage_table_post, cat_coverage_table_cumulative


rule all:
	input:
		expand('downsampled/downsampled_{fname}_{cov}x.cov.gz', fname=merged_reps, cov=covgs1),
		expand('downsampled/downsampled_{fname}_{cov}x.cov.gz', fname=merged_samples, cov=covgs2),
		'total_coverage/coverage_table.txt',
		#'cumulative_coverage/cumulative_coverage_table.txt',
		#'total_coverage/coverage_table_post.txt'

rule determine_coverage:
	input:
		os.path.join(bed_dir, '{fname}.cov.gz')
	output:
		'total_coverage/{fname}_total_coverage.txt'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < file_coverage.R
		"""

rule cat_coverage_table:
	input:
		expand('total_coverage/{fname}_total_coverage.txt', fname=merged_reps + merged_samples)
	output:
		'total_coverage/coverage_table.txt',

	shell:
		"""
		cat {input} > {output[0]}
		"""

def get_coverage1(wildcards, cov_file):
	covs = []
	sample_coverage = None
	for line in open(cov_file, 'r'):
		prep = line.split()[0]
		sample = line.split()[1]
		coverage = line.split()[2]
		if sample in wildcards.sample:
			covs.append(int(coverage))
			if prep in wildcards.prep:
				sample_coverage = coverage

	min_cov = min(covs)
	fraction = round(float(min_cov)/float(sample_coverage), 5)
	return(fraction)

def get_coverage2(wildcards, cov_file):
	sample_coverage = None
	for line in open(cov_file, 'r'):
		prep = line.split()[0]
		sample = line.split()[1]
		coverage = line.split()[3]
		if sample in wildcards.sample and prep in wildcards.prep:
				sample_coverage = coverage
				break
	fraction = round(float(wildcards.cov)/float(sample_coverage), 5)
	return(fraction)

rule downsample_bedgraph:
	input:
		os.path.join(bed_dir, 'repsCombined_{prep}_{sample}.cov.gz'),
		'total_coverage/coverage_table.txt'

	params:
		sample_fraction = lambda wildcards, input: get_coverage2(wildcards, input[1])

	output:
		#'downsampled/downsampled_{prep}_{sample}_' + params.min_coverage + '.cov.gz'
		'downsampled_{cov}x/downsampled_{prep}_{sample}_{cov}x.cov.gz'

	run:
		if params.sample_fraction < 1:
			shell('zcat {input[0]} | python3 downsample_methylKit.py --fraction {params.sample_fraction} --bedGraph | gzip > {output}')
		else:
			#shell('cp {input[0]} {output}')
			shell('touch {output}')



rule determine_coverage_post:
	input:
		'downsampled_{cov}x/downsampled_{prep}_{sample}_{cov}x.cov.gz'
	output:
		'total_coverage/{prep}_{sample}_{cov}_total_coverage.txt'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < file_coverage.R
		"""

rule cat_coverage_table_post:
	input:
		expand('total_coverage/{prep}_{sample}_{cov}_total_coverage.txt', prep=preps, sample=samples, cov=covgs)
	output:
		'total_coverage/coverage_table_post.txt',
	params:
		preps = ','.join(preps),
		samples = ','.join(samples)

	shell:
		"""
		cat {input} > {output[0]}
		"""

rule cumulative_coverage:
	input:
		'downsampled_{cov}x/downsampled_{prep}_{sample}_{cov}x.cov.gz'
	output:
		'cumulative_coverage/{prep}_{sample}_{cov}_cumulative_coverage.txt'
	shell:
		"""
		module load R
		R --no-restore --no-save --args {input} {output} < coverage_table.R
		"""

rule cat_coverage_table_cumulative:
	input:
		expand('cumulative_coverage/{prep}_{sample}_{cov}_cumulative_coverage.txt', prep=preps, sample=samples, cov=covgs)
	output:
		'cumulative_coverage/cumulative_coverage_table.txt',
	params:
		preps = ','.join(preps),
		samples = ','.join(samples)

	shell:
		"""
		cat {input} > {output[0]}
		"""
