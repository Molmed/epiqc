#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --time=168:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=EpiQC_bsseeker2
#SBATCH --output=logs/EpiQC_bsseeker2.SAMPLE.%A_%a.log

echo [`date`] "Started job"

### metadata and refs
outdir=seqc2/epiqc/bsseeker2/outs
refdir=seqc2/epiqc/bsseeker2/grch38_core_and_bs_controls.fa_bowtie2
ref=seqc2/epiqc/grch38_core_bs_controls/grch38_core_and_bs_controls.fa

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
if [ ! -f ${outdir}/${sample}.sorted.bam ] && [ ! -f ${outdir}/${sample}.sorted.md.bam ]; then
  rsync -avL $R1path $TMPDIR
  rsync -avL $R2path $TMPDIR
  
  ## adapters already trimmed
  R1trim=$R1
  R2trim=$R2
  
  ## alignment
  conda activate bsseeker2  # a python2 enviroment
  bs_seeker2-align.py --aligner=bowtie2 -1 $R1trim -2 $R2trim -g $ref -o ${sample}.bam --bt-p 16
  sambamba sort -m 12G -t 8 ${sample}.bam
  rsync -av ${sample}.sorted.bam* $outdir
fi

### deduplicate BAM
if [ ! -f ${outdir}/${sample}.sorted.md.bam ]; then
  rsync -av ${outdir}/${sample}.sorted.bam* $TMPDIR
  java -jar -Xmx12G picard.jar MarkDuplicates \
      I=${sample}.sorted.bam \
      O=${sample}.sorted.md.bam \
      CREATE_INDEX=true \
      M=dedup_metrics_${sample}.txt
  rm ${sample}.sorted.bam
  rsync -av ${sample}.sorted.md.bam* $outdir
fi

### call methylation
if [ ! -f ${outdir}/${sample}.CGmap.gz ]; then
  conda activate bsseeker2
  if [ ! -f ${sample}.sorted.md.bam ] ; then rsync -av ${outdir}/${sample}.sorted.md.bam* $TMPDIR; fi
  bs_seeker2-call_methylation.py -i ${sample}.sorted.md.bam -d $refdir --sorted -o ${sample}
  rm $R1trim $R2trim
  rsync -av ${sample}*.gz $outdir
  rsync -av ${sample}*log $outdir
fi

echo [`date`] "Finished job"
