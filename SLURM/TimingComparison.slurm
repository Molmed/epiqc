#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=32G
#SBATCH --job-name=EpiQC_timecomp
#SBATCH --output=logs/EpiQC_timecomp.%A_%a.log

source ~/.bashrc
echo [`date`] "Started job"

### metadata and refs
rsync -av seqc2/epiqc/grch38_core_bs_controls/grch38_core_and_bs_controls* $TMPDIR
ref=grch38_core_and_bs_controls.fa
samplesheet=$1 # a file with absolute path to FASTQs (R1), one per line (in this case, subsampled to 1M reads using seqtk sample -s1234)
algo=$2  # can be bwmeth, bismark, bitmapperBS, or bsseeker2
R1path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
R2path=${R1path/.R1/.R2}
rsync -av $R1path $TMPDIR
rsync -av $R2path $TMPDIR
R1=$(basename $R1path)
R2=$(basename $R2path)
sample=$(echo $R1 | cut -d'.' -f1)
echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "This job is for algo: " $algo
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### iterate ten times (parse the logs files to get timings)
i=0
while (( $i < 10 )); do
  if [[ $algo == "bwameth" ]]; then
    conda activate bwameth
    time bwameth.py -t 4 --reference $ref $R1 $R2 > asdf.sam 
    sambamba view -t 2 -S -f bam -o "${sample}.bam" asdf.sam
    sambamba sort ${sample}.bam
    time MethylDackel extract -@4 $ref ${sample}.sorted.bam
    ((i+=1))
  elif [[ $algo == "bismark" ]]; then
    conda activate bwameth
    time bismark --parallel 2 --genome seqc2/epiqc/grch38_core_bs_controls -1 $R1 -2 $R2
    time bismark_methylation_extractor --gzip --multicore 2 --bedGraph --zero_based ${sample}.R1.1M_bismark_bt2_pe.bam
    ((i+=1))
  elif [[ $algo == "bitmapperBS" ]]; then
    conda activate
    time bitmapperBS --search seqc2/epiqc/bitmapperbs/ref_genome --seq1 $R1 --seq2 $R2 --bam --mapstats ${sample}.mapstats -o ${sample}.bam
    sambamba sort ${sample}.bam
    time MethylDackel extract -@4 $ref ${sample}.sorted.bam
    ((i+=1))
  elif [[ $algo == "bsseeker2" ]]; then
    conda activate bsseeker2
    time /bs_seeker2-align.py --aligner=bowtie2 -1 $R1 -2 $R2 -g $ref -o ${sample}.bam --bt-p 4
    sambamba sort ${sample}.bam
    time bs_seeker2-call_methylation.py -i ${sample}.sorted.bam -d seqc2/epiqc/bsseeker2/grch38_core_and_bs_controls.fa_bowtie2 --sorted -o ${sample}
    ((i+=1))
  else
    echo "something went wrong with $algo"
  fi
  conda deactivate
done
echo [`date`] "Finished job"
