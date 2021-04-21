#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00
#SBATCH --mem=32G
#SBATCH --job-name=EpiQC_bwameth
#SBATCH --output=logs/EpiQC_bwameth.%A_%a.log

echo [`date`] "Started job"

### metadata and refs
outdir=seqc2/epiqc/bwameth/outs
rsync -av seqc2/epiqc/grch38_core_bs_controls/grch38_core_and_bs_controls* $TMPDIR
ref=grch38_core_and_bs_controls.fa
samplesheet=$1  # a file with absolute path to one fastq (R1) per line

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
if [ ! -f ${outdir}/${sample}.bam ] && [ ! -f ${outdir}/${sample}.sorted.bam ] && [ ! -f ${outdir}/${sample}.sorted.md.bam ]; then
  rsync -avL $R1path $TMPDIR
  rsync -avL $R2path $TMPDIR
  sampleRG=$(echo $sample | cut -d'-' -f1)
  readgroup=@RG\\tID:${sample}\\tSM:${sampleRG}\\tPL:ILLUMINA
  
  ## adapters already trimmed
  R1trim=$R1
  R2trim=$R2
  
  ## run bwameth
  bwameth.py -t 8 --read-group $readgroup --reference $ref $R1trim $R2trim \
  | sambamba view -t 2 -S -f bam -o "${sample}.bam" /dev/stdin
  rsync -av ${sample}.bam* $outdir
fi

### if needed, create sorted (and indexed) BAM
if [ ! -f ${outdir}/${sample}.sorted.bam ] && [ ! -f ${outdir}/${sample}.sorted.md.bam ]; then
  rsync -av ${outdir}/${sample}.bam $TMPDIR
  sambamba sort -m 12G -t 1 ${sample}.bam
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
if [ ! -f ${outdir}/${sample}.5mC.bedGraph ]; then
  rsync -av ${outdir}/${sample}.sorted.md.bam* $TMPDIR
  echo "[`date`] Calling methylation" && MethylDackel extract -@4 $ref ${sample}.sorted.md.bam
  mv ${sample}.sorted.md_CpG.bedGraph ${sample}.5mC.bedGraph
  rsync -av ${sample}.5mC.bedGraph $outdir
fi

### merge contexts
if [ ! -f ${outdir}/${sample}.5mC.merged.bedGraph* ]; then
  rsync -av ${outdir}/${sample}.5mC.bedGraph $TMPDIR
  echo "[`date`] Merging strands" && MethylDackel mergeContext $ref ${sample}.5mC.bedGraph > ${sample}.5mC.merged.bedGraph
  echo "[`date`] Gzipping" && gzip ${sample}.5mC.merged.bedGraph
  rsync -av ${sample}.5mC.merged.bedGraph.gz $outdir
fi

### create methylKit format
if [ ! -f ${outdir}/${sample}.5mC.merged.methylKit* ]; then
  rsync -av ${outdir}/${sample}.5mC.merged.bedGraph.gz $TMPDIR
  echo "[`date`] Converting to methylKit" && methTrans_biscov2methylKit.py ${TMPDIR}/${sample}.5mC.merged.bedGraph.gz
  rsync -av ${sample}.5mC.merged.methylKit.gz $outdir
fi

echo [`date`] "Finished job"
