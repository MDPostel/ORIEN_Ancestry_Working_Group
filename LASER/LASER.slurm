#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=220G
#SBATCH --cpus-per-task=11
#SBATCH --ntasks=20
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=laser

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
module load usc bzip2 bcftools htslib vcftools
t=$(echo 4000)
ulimit -n 4001

HGDP=$(echo HGDP.vcf.gz)
ORIEN=$(echo ORIEN.vcf.gz)
GP1000=$(echo 1000GP.vcf.gz)
TCGA=$(echo TCGA.vcf.gz)

bcftools merge $ORIEN $TCGA --threads 4000 -Oz -o ORIEN_TCGA.vcf.gz
tabix -f ORIEN_TCGA.vcf.gz

#in trace.conf:
#STUDY_FILE      ORIEN_TCGA_laser_samples.geno
#GENO_FILE       1000GP_HGDP_laser_refs.geno
#COORD_FILE      /project/davidwcr_788/ORIEN_Intermember_Project/ancestry_tests/LASER/epyc-64/laser.RefPC.coord

/home1/mpostel/bin/LASER-2.04/vcf2geno/vcf2geno --inVcf ORIEN_TCGA.vcf.gz --out ORIEN_TCGA_laser_samples
/home1/mpostel/bin/LASER-2.04/vcf2geno/vcf2geno --inVcf HGDP_1000GP.vcf.gz --out 1000GP_HGDP_laser_refs

cd /project/davidwcr_788/ORIEN_Intermember_Project/ancestry_tests/LASER/epyc-64

#create coord file so program doesn't have to rerun it each time
geno=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/1000GP_HGDP_laser_refs.geno)
/home1/mpostel/bin/LASER-2.04/laser -g $geno -pca 1 -k 6

#Run in parallel in one script or create multiple scripts that all reference the coord file made above
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 1 -y 1250 -o results.1-1250 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 1251 -y 2500 -o results.1251-2500 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 2501 -y 3750 -o results.2501-3750 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 3751 -y 5000 -o results.3751-5000 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 5001 -y 6250 -o results.5001-6250 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 6251 -y 7500 -o results.6251-7500 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 7501 -y 8750 -o results.7501-8750 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 8751 -y 10000 -o results.8751-10000 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 10001 -y 11250 -o results.10001-11250 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 11251 -y 12500 -o results.11251-12500 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 12501 -y 13750 -o results.12501-13750 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 13751 -y 15000 -o results.13751-15000 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 15001 -y 16250 -o results.15001-16250 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 16251 -y 17500 -o results.16251-17500 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 17501 -y 18750 -o results.17501-18750 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 18751 -y 20000 -o results.18751-20000 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 20001 -y 21250 -o results.20001-21250 &
/home1/mpostel/bin/LASER-2.04/trace -p trace.conf -x 21251 -y 22500 -o results.21251-21608 &
