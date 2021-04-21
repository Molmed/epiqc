#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=EpiQC_adaptertrim
#SBATCH --output=logs/EpiQC_adaptertrim.%A_%a.log

source ~/.bashrc
conda activate bwameth
echo [`date`] "Started job"

### metadata and refs
outdir=/athena/masonlab/scratch/projects/seqc2/epiqc/revision/data_adaptertrimmed/outs
samplesheet=$1
R1path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
R2path=${R1path/.R1/.R2}
R1=$(basename $R1path)
R2=$(basename $R2path)
sample=$(echo $R1 | cut -d'.' -f1)
echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### initial alignment
if [ ! -f ${outdir}/${sample}.R1.trimmed.fastq.gz ] || [ ! -f ${outdir}/${sample}.R2.trimmed.fastq.gz ]; then
  rsync -avL $R1path $TMPDIR
  rsync -avL $R2path $TMPDIR
  
  ## trim adapters (all non-MethylSeq)
  if [[ $sample != "MethylSeq"* ]]; then
  fastp --in1 $R1 --in2 $R2 \
    --out1 ${R1/.fastq.gz/.trimmed.fastq.gz} --out2 ${R2/.fastq.gz/.trimmed.fastq.gz} \
    -l 2 -Q -z 4 --trim_poly_g --overrepresentation_analysis

  ## MethylSeq needs Adaptase hard trimmed
  else
    fastp --in1 $R1 --in2 $R2 \
    --out1 ${R1/.fastq.gz/.trimmed.fastq.gz} --out2 ${R2/.fastq.gz/.trimmed.fastq.gz} \
    --trim_front1 10 --trim_tail1 10 --trim_front2 10 \
    -l 2 -Q -z 4 --trim_poly_g --overrepresentation_analysis
  fi

  # export
  R1trim=${R1/.fastq.gz/.trimmed.fastq.gz}
  R2trim=${R2/.fastq.gz/.trimmed.fastq.gz}
  rsync -av $R1trim $outdir
  rsync -av $R2trim $outdir
fi

echo [`date`] "Ended job"
