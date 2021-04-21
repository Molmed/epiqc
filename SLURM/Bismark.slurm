#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=EpiQC_bismark 
#SBATCH --output=logs/EpiQC_bismark.%A_%a.log

echo [`date`] 'Started job'

### metadata and refs
outdir=seqc2/epiqc/bismark/outs
samplesheet=$1  # a file with the absolute path to one fastq (R1) per line
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
if [ ! -f ${outdir}/${sample}.R1.trimmed_bismark_bt2_pe.bam ] && [ ! -f ${outdir}/${sample}.R1.trimmed_bismark_bt2_pe.sorted.bam ] && [ ! -f ${outdir}/${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam ]; then
  rsync -avL $R1path $TMPDIR
  rsync -avL $R2path $TMPDIR
  sampleSM=$(echo $sample | cut -d'-' -f1)
  
  ## adapters already trimmed
  R1trim=$R1
  R2trim=$R2
    
  ## run BISMARK
  bismark --parallel 2 --nucleotide_coverage \
    --rg_tag --rg_id $sample --rg_sample $sampleSM \
    --genome seqc2/epiqc/grch38_core_bs_controls \
    -1 $R1trim \
    -2 $R2trim \
    --temp_dir /athena/masonlab/scratch/users/jfoox/tmp/${sample}
  rsync -av ${sample}*bismark* $outdir

  ## sort BAM
  sambamba sort -n -m 12G -t 4 -p ${sample}.R1.trimmed_bismark_bt2_pe.bam
  
  ## export (once above for safety)
  rsync -av ${sample}*bismark* $outdir
fi

### deduplicate BAM
if [ ! -f ${outdir}/${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam ]; then
  if [ ! -f ${sample}.R1.trimmed_bismark_bt2_pe.sorted.bam ]; then rsync -av ${outdir}/${sample}.R1.trimmed_bismark_bt2_pe.sorted.bam* $TMPDIR; fi
  deduplicate_bismark --bam ${sample}.R1.trimmed_bismark_bt2_pe.sorted.bam

  ## export/index and export
  rsync -av ${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam $outdir
  sambamba index -t 4 ${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam
  rsync -av ${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam.bai $outdir
  rsync -av *.txt $outdir
fi

### methylation extraction
if [ ! -f ${outdir}/${sample}*bismark.cov.gz ]; then
  if [ ! -f ${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam ]; then rsync -av ${outdir}/${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam* $TMPDIR; fi
  bismark_methylation_extractor --gzip --multicore 2 --bedGraph --zero_based ${sample}.R1.trimmed_bismark_bt2_pe.sorted.deduplicated.bam
  rsync -av ${sample}*.bismark.cov.gz $outdir
fi

echo [`date`] "Finished job"
