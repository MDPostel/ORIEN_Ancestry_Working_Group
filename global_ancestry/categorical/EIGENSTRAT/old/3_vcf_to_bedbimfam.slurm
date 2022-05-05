#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=plink

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen
gunzip concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP.vcf.gz
#vcftools --vcf concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP.vcf --plink --out concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP
plink --vcf concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP.vcf --recode --out concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP_plink --double-id
