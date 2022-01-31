#!/usr/bin/bash
#SBATCH --partition=epyc-64
##SBATCH --exclusive
#SBATCH --mem=240G
##SBATCH --cpus-per-task=16
##SBATCH --ntasks=64
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=admixture

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
isec=$(echo HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz)

~/bin/plink-1.90/plink --vcf $isec --double-id --make-bed --out isec
bed=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/isec.bed)
K=$(echo 8)
threads=$(echo 4000)
cd /project/davidwcr_788/ORIEN_Intermember_Project/ancestry_tests/admixture 
/home1/mpostel/bin/admixture/admixture_linux-1.3.0/admixture ${bed} ${K} -j${threads}

