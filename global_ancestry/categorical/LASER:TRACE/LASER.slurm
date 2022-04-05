#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=laser

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
module load usc bzip2 bcftools htslib vcftools
merge4=$(echo merge.recode.vcf.gz)
t=$(echo 4000)
ulimit -n 4001

#from AIMs script:
        #bcftools view $merge4 --force-samples -S HGDP_sample_names.txt --threads $t > HGDP_postmerge.vcf
        #bcftools view $merge4 --force-samples -S ORIEN_sample_names.txt --threads $t > ORIEN_postmerge.vcf
        #bcftools view $merge4 --force-samples -S 1000GP_sample_names.txt --threads $t > 1000GP_postmerge.vcf
        #bcftools view $merge4 --force-samples -S TCGA_sample_names.txt --threads $t > TCGA_postmerge.vcf
        #bgzip HGDP_postmerge.vcf; tabix -f HGDP_postmerge.vcf.gz
        #bgzip ORIEN_postmerge.vcf; tabix -f ORIEN_postmerge.vcf.gz
        #bgzip 1000GP_postmerge.vcf; tabix -f 1000GP_postmerge.vcf.gz
        #bgzip TCGA_postmerge.vcf; tabix -f TCGA_postmerge.vcf.gz
        #cat HGDP_sample_names.txt > HGDP_1000GP_samples.txt
        #cat 1000GP_sample_names.txt >> HGDP_1000GP_samples.txt
        #bcftools view $merge4 --force-samples -S HGDP_1000GP_samples.txt --threads $t -Oz -o HGDP_1000GP.vcf.gz
        #tabix -f HGDP_1000GP.vcf.gz
HGDP=$(echo HGDP_postmerge.vcf.gz)
ORIEN=$(echo ORIEN_postmerge.vcf.gz)
GP1000=$(echo 1000GP_postmerge.vcf.gz)
TCGA=$(echo TCGA_postmerge.vcf.gz)

cat ORIEN_sample_names.txt > ORIEN_TCGA_samples.txt
cat TCGA_sample_names.txt >> ORIEN_TCGA_samples.txt
bcftools view $merge4 --force-samples -S ORIEN_TCGA_samples.txt --threads $t -Oz -o ORIEN_TCGA.vcf.gz
tabix -f ORIEN_TCGA.vcf.gz

#in trace.conf:
#STUDY_FILE      ORIEN_TCGA_laser_samples.geno
#GENO_FILE       1000GP_HGDP_laser_refs.geno

/home1/mpostel/bin/LASER-2.04/vcf2geno/vcf2geno --inVcf ORIEN_TCGA.vcf.gz --out ORIEN_TCGA_laser_samples
/home1/mpostel/bin/LASER-2.04/vcf2geno/vcf2geno --inVcf HGDP_1000GP.vcf.gz --out 1000GP_HGDP_laser_refs

/home1/mpostel/bin/LASER-2.04/trace -p trace.conf
