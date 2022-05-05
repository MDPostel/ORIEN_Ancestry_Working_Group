#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=400GB
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=lamp


#echo "populations=3" > config.txt
#echo "genofile=ORIEN_USC_samples_merged3.geno" >> config.txt
#echo "posfile=/home1/mpostel/bin/LAMP2.5/test/chr1.pos" >> config.txt

lamp config.txt


