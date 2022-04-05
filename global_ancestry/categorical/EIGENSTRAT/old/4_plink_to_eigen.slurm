#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=convertf

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen
convertf -p parfile
