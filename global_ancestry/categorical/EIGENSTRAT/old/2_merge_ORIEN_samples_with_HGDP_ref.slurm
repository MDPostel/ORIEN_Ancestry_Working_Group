
#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=vcfmerge

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen

ulimit -n 3000

echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_samples_merged.vcf" > gVCF_USCORIEN_HGDP_list.txt; echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/HGDP_samples_merged.vcf" >> gVCF_USCORIEN_HGDP_list.txt

bcftools merge --file-list gVCF_USCORIEN_HGDP_list.txt -g /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_num_without_chr.fa -O v -o concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP.vcf

bgzip concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP.vcf
tabix concat_WES_snps_for_ancestry_gVCF_ORIEN_HGDP.vcf.gz
