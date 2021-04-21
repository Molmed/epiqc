#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=32G
#SBATCH --job-name=EpiQC_algo_mBias
#SBATCH --output=logs/EpiQC_algo_mBias.%A_%a.log

echo [`date`] Started job

### metadata and refs
outdir=seqc2/epiqc/mbias/outs
samplesheet=$1  # a file with absolute path to BAMs, one per line
bampath=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $bampath $TMPDIR
bam=$(basename $bampath)
sample=$(echo $bam | cut -d'.' -f1)
algo=$(echo $bampath | cut -d'/' -f9) # e.g. bismark
rsync -av seqc2/epiqc/grch38_core_bs_controls/grch38_core_and_bs_controls.fa* $TMPDIR
ref=grch38_core_and_bs_controls.fa
echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "/////////////////////////////////////////////////"

### run
cd $TMPDIR

### adjust bismark from position sorted to coordinate sorted
if [[ $algo == "bismark" ]]; then
  sambamba sort -t 4 $bam
  MethylDackel mbias -@ 4 --txt $ref ${bam/.bam/.sorted.bam} ${algo}_${sample} > mBias_table_${algo}_${sample}.txt 2> mBias_suggestion_${algo}_${sample}.txt
else
  MethylDackel mbias -@ 4 --txt $ref $bam ${algo}_${sample} > mBias_table_${algo}_${sample}.txt 2> mBias_suggestion_${algo}_${sample}.txt
fi

### export
rsync -av *.txt $outdir

echo [`date`] Finished job

