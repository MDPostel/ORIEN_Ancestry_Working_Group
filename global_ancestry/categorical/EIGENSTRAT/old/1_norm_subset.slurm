#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=subset

echo "ORIEN Samples first..."
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_gVCFs
for file in `ls *vcf.gz | grep -v FIX`
    do
        echo "working on" $file
                name=`echo $file | cut -d. -f1`
        bcftools norm -f /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa $file > ${name}_norm.vcf
        bgzip ${name}_norm.vcf
        tabix ${name}_norm.vcf.gz
        bcftools annotate --rename-chrs /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/chr_name_conv.txt ${name}_norm.vcf.gz -Oz -o ${name}_norm_chrFIX.vcf.gz
        gunzip ${name}_norm_chrFIX.vcf.gz
        bgzip ${name}_norm_chrFIX.vcf
        tabix ${name}_norm_chrFIX.vcf.gz
        bcftools view --types snps -R /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/170_SNPs_list_altlociremoved_ordered.txt ${name}_norm_chrFIX.vcf.gz > /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/${name}_bedsplit.vcf
        bgzip /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/${name}_bedsplit.vcf
        tabix /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/${name}_bedsplit.vcf.gz
done

echo "merging ORIEN sample SNP gVCFs"
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/
ls /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/*_bedsplit.vcf.gz > samples_list.txt
bcftools merge --file-list samples_list.txt -g /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_num_without_chr.fa -O v -o ORIEN_USC_samples_merged.vcf

echo "Moving on to HGDP..."
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/HGDP_ref_files
for file in `ls *wgs*vcf.gz`
    do
        echo "working on" $file
                name=`echo $file | cut -d. -f1`
        bcftools norm -f /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa $file > ${name}_norm.vcf
        bgzip ${name}_norm.vcf
        tabix ${name}_norm.vcf.gz
        bcftools annotate --rename-chrs /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/chr_name_conv.txt ${name}_norm.vcf.gz -Oz -o ${name}_norm_chrFIX.vcf.gz
        gunzip ${name}_norm_chrFIX.vcf.gz
        bgzip ${name}_norm_chrFIX.vcf
        tabix ${name}_norm_chrFIX.vcf.gz
        bcftools view --types snps -R /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/170_SNPs_list_altlociremoved_ordered.txt ${name}_norm_chrFIX.vcf.gz > /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/${name}_HGDPref_bedsplit.vcf
        bgzip /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/${name}_HGDPref_bedsplit.vcf
        tabix /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/${name}_HGDPref_bedsplit.vcf.gz
done

echo "merging HGDP SNP gVCFs"
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/
ls /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/ORIEN_USC_eigen/bed_split_gVCFs/*_HGDPref_bedsplit.vcf.gz > ref_list.txt
bcftools merge --file-list ref_list.txt  -g /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_num_without_chr.fa -O v -o HGDP_samples_merged.vcf
