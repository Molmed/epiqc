#!/bin/bash
#SBATCH -A snic2017-7-344
##SBATCH -t 12:00:00
#SBATCH -t 01:00:00
#SBATCH -p core -n 1
#SBATCH -J snakemake

module load bioinfo-tools snakemake

snakemake -j 100 --cluster-config downsample.json --cluster "sbatch -A {cluster.A} -t {cluster.t} -p {cluster.p} -n {cluster.n}" -s downsample.smk.py
