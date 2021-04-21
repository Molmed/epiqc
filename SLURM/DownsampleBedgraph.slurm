#!/bin/bash -l
#SBATCH --time=168:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=EpiQC_downsampleBedgraphs
#SBATCH --output=logs/EpiQC_downsampleBedgraphs.SAMPLE.%A_%a.log

echo [`date`] Started job

# set up
samplesheet=$1  # a file with the path to one bedGraph on each line
bedgraphdir=$2
outdir=/athena/masonlab/scratch/projects/seqc2/epiqc/revision/bedgraphs/outs
sample=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av ${bedgraphdir}/${sample}.bedGraph.gz $TMPDIR
rsync -av ${bedgraphdir}/downsample_NEB.py $TMPDIR
cd $TMPDIR

cov=$(cat ${bedgraphdir}/${sample}.meancov.txt)  # calculated per sample with: cat ${sample}.bedGraph | awk '{ sum += $5 + $6 } END { print sum / NR }' > ${sample}.meancov.txt
for dsgoal in 5 10 15 20 25 30 40; do 
  fraction=$(echo "scale=3; $dsgoal / $cov" | bc)
  bg=${sample}.bedGraph.gz
  
  roundedcov=$(echo $cov | awk '{print int($1+0.5)}')  # e.g. 24.6 --> 25
  if [ "$roundedcov" -lt "$dsgoal" ]; then
    echo "DS target (${dsgoal}) larger than average cov for " $bg
  else
    fraction=$(echo "scale=3; $dsgoal / $cov" | bc)
    zcat $bg | grep -v track | python downsample_NEB.py --fraction $fraction --bedGraph > ${sample}.ds${dsgoal}.bedGraph
    rsync -av ${sample}.ds${dsgoal}.bedGraph $outdir
  fi
done


echo [`date`] Ended job
