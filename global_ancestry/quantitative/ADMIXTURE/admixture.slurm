#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=240G
#SBATCH --cpus-per-task=64
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=admixsupervised

threads=$SLURM_CPUS_PER_TASK
${admixture} ${refs} ${K} -j${threads} --cv=10

cp HGDP_1000GP.${K}.P ORIEN_TCGA.${K}.P.in
${admixture} -P ${samples} ${K} -j${threads} --cv=10
