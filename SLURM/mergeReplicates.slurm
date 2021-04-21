#!/bin/bash -l
#SBATCH --time=168:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=EpiQC_mergeBedgraphs
#SBATCH --output=logs/EpiQC_mergeBedgraphs.%A_%a.log

echo [`date`] Started job

# set up
samplesheet=$1  # a file with absolute path to bedGraphs, one per line
bedgraphdir=$2
outdir=seqc2/epiqc/revision/bedgraphs/outs
sample=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
cd $TMPDIR

if [ -f ${bedgraphdir}/${sample} ]

echo [`date`] Merging ...
### TrueMethyl is labeled slightly differently
if [[ $sample == "TrueMethyl"* ]]; then
  genome=$(echo $sample | cut -d_ -f2)
  rsync -av ${bedgraphdir}/TrueMethyl*${genome}*.5mC.merged.bedGraph.gz $TMPDIR
  combine_methylation_outputs.py $sample TrueMethyl*.bedGraph.gz

### all other cases
else
  rsync -av ${bedgraphdir}/${sample}*.5mC.merged.bedGraph.gz $TMPDIR
  combine_methylation_outputs.py $sample ${sample}*.bedGraph.gz
fi

# get mean
zcat ${sample}.bedGraph.gz | awk ' { sum += $5+$6 } END { print sum / NR } ' > ${sample}.meancov.txt
# get median
zcat ${sample}.bedGraph.gz | awk ' { print $5+$6 } ' | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' > ${sample}.mediancov.txt

# transfer
rsync -av ${sample}.bedGraph.gz $outdir
rsync -av ${sample}.*cov.txt $outdir
