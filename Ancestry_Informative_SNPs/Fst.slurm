#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=850G
#SBATCH --time=07:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=Fst

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
vcf=$(echo "HGDP_1000GP_ORIEN_minusTCGA_straight2plink_bed2vcf_addchr.vcf.gz")
module load usc bzip2 bcftools htslib vcftools
plink=$(echo "/home1/mpostel/bin/plink-1.90/plink")

#Download lists of HGDP and 1000GP individuals from: https://www.internationalgenome.org/data-portal/sample
#rename to:
#igsr_samples_HGDP.tsv
#igsr_samples_1000GP.tsv

#awk -F '\t' -v OFS='\t' '{print $1,$7}' igsr_samples_HGDP.tsv > SampleID_Pop_HGDP.tsv
#awk -F '\t' -v OFS='\t' '{print $1,$6}' igsr_samples_1000GP.tsv > SampleID_Pop_1000GP.tsv

#awk -F '\t' '{print $2}' SampleID_Pop_HGDP.tsv | sort -u | grep -v Superpopulation > HGDP_pops.txt
#awk -F '\t' '{print $2}' SampleID_Pop_1000GP.tsv | sort -u | grep -v Superpopulation > 1000GP_pops.txt
#sed -i 's/Africa (HGDP)/AFR/g' HGDP_pops.txt; sed -i 's/America (HGDP)/AMR/g' HGDP_pops.txt; sed -i 's/Central South Asia (HGDP)/SAS/g' HGDP_pops.txt; sed -i 's/East Asia (HGDP)/EAS/g' HGDP_pops.txt; sed -i 's/Europe (HGDP)/EUR/g' HGDP_pops.txt; sed -i 's/Middle East (HGDP)/MDE/g' HGDP_pops.txt; sed -i 's/Oceania (SGDP)\,Oceania (HGDP)/OCE/g' HGDP_pops.txt; sed -i 's/Oceania (HGDP)/OCE/g' HGDP_pops.txt
#sed -i 's/Africa (HGDP)/AFR/g' SampleID_Pop_HGDP.tsv; sed -i 's/America (HGDP)/AMR/g' SampleID_Pop_HGDP.tsv; sed -i 's/Central South Asia (HGDP)/SAS/g' SampleID_Pop_HGDP.tsv; sed -i 's/East Asia (HGDP)/EAS/g' SampleID_Pop_HGDP.tsv; sed -i 's/Europe (HGDP)/EUR/g' SampleID_Pop_HGDP.tsv; sed -i 's/Middle East (HGDP)/MDE/g' SampleID_Pop_HGDP.tsv; sed -i 's/Oceania (SGDP)\,Oceania (HGDP)/OCE/g' SampleID_Pop_HGDP.tsv; sed -i 's/Oceania (HGDP)/OCE/g' SampleID_Pop_HGDP.tsv

#grep AFR SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > AFR.txt; grep AFR SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> AFR.txt
#grep AMR SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > AMR.txt; grep AMR SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> AMR.txt
#grep SAS SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > SAS.txt; grep SAS SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> SAS.txt
#grep EAS SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > EAS.txt; grep EAS SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> EAS.txt
#grep EUR SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > EUR.txt; grep EUR SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> EUR.txt
#grep OCE SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > OCE.txt; grep OCE SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> OCE.txt
#grep MDE SampleID_Pop_1000GP.tsv | grep -oe HG[0-9]* -oe NA[0-9]* > MDE.txt; grep MDE SampleID_Pop_HGDP.tsv | grep -o HGDP[0-9]* >> MDE.txt

#step above makes 5 super pops + Oceania + Middle East:
#AFR
#AMR
#SAS
#EAS
#EUR
#OCE = oceania
#MDE = middle east

#for all txt files -- you may have to do ```for i in ${cat file.txt}; do echo ${i}_${i}; done``` because plink doubled ID names to be both indiv and family

#cat AMR_indivfam.txt > notAFR_indivfam.txt; cat EAS_indivfam.txt >> notAFR_indivfam.txt; cat EUR_indivfam.txt >> notAFR_indivfam.txt; cat MDE_indivfam.txt >> notAFR_indivfam.txt; cat OCE_indivfam.txt >> notAFR_indivfam.txt; cat SAS_indivfam.txt >> notAFR_indivfam.txt
#cat AFR_indivfam.txt > notAMR_indivfam.txt; cat EAS_indivfam.txt >> notAMR_indivfam.txt; cat EUR_indivfam.txt >> notAMR_indivfam.txt; cat MDE_indivfam.txt >> notAMR_indivfam.txt; cat OCE_indivfam.txt >> notAMR_indivfam.txt; cat SAS_indivfam.txt >> notAMR_indivfam.txt
#cat AFR_indivfam.txt > notEAS_indivfam.txt; cat AMR_indivfam.txt >> notEAS_indivfam.txt; cat EUR_indivfam.txt >> notEAS_indivfam.txt; cat MDE_indivfam.txt >> notEAS_indivfam.txt; cat OCE_indivfam.txt >> notEAS_indivfam.txt; cat SAS_indivfam.txt >> notEAS_indivfam.txt
#cat AFR_indivfam.txt > notEUR_indivfam.txt; cat AMR_indivfam.txt >> notEUR_indivfam.txt; cat EAS_indivfam.txt >> notEUR_indivfam.txt; cat MDE_indivfam.txt >> notEUR_indivfam.txt; cat OCE_indivfam.txt >> notEUR_indivfam.txt; cat SAS_indivfam.txt >> notEUR_indivfam.txt
#cat AFR_indivfam.txt > notMDE_indivfam.txt; cat AMR_indivfam.txt >> notMDE_indivfam.txt; cat EAS_indivfam.txt >> notMDE_indivfam.txt; cat EUR_indivfam.txt >> notMDE_indivfam.txt; cat OCE_indivfam.txt >> notMDE_indivfam.txt; cat SAS_indivfam.txt >> notMDE_indivfam.txt
#cat AFR_indivfam.txt > notOCE_indivfam.txt; cat AMR_indivfam.txt >> notOCE_indivfam.txt; cat EAS_indivfam.txt >> notOCE_indivfam.txt; cat EUR_indivfam.txt >> notOCE_indivfam.txt; cat MDE_indivfam.txt >> notOCE_indivfam.txt; cat SAS_indivfam.txt >> notOCE_indivfam.txt
#cat AFR_indivfam.txt > notSAS_indivfam.txt; cat AMR_indivfam.txt >> notSAS_indivfam.txt; cat EAS_indivfam.txt >> notSAS_indivfam.txt; cat EUR_indivfam.txt >> notSAS_indivfam.txt; cat OCE_indivfam.txt >> notSAS_indivfam.txt; cat MDE_indivfam.txt >> notSAS_indivfam.txt

#--------

bcftools view $vcf -S HGDP_and_1000GP_sample_names_indivfam.txt --threads 4000 -Oz -o HGDP_1000GP_onlyrefpops.vcf.gz; tabix -f HGDP_1000GP_onlyrefpops.vcf.gz

vcf2=$(echo "HGDP_1000GP_onlyrefpops.vcf.gz")

vcftools --gzvcf HGDP_1000GP.vcf.gz --weir-fst-pop AFR.txt --weir-fst-pop AMR.txt --weir-fst-pop SAS.txt --weir-fst-pop EAS.txt --weir-fst-pop EUR.txt --weir-fst-pop OCE.txt --weir-fst-pop MDE.txt --out overall_fsts

vcftools --gzvcf HGDP_1000GP.vcf.gz --weir-fst-pop AFR.txt --weir-fst-pop AMR.txt --weir-fst-pop SAS.txt --weir-fst-pop EAS.txt --weir-fst-pop EUR.txt --out overall_fsts_5trunks_no_OCE_no_MDE

vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop notAFR_indivfam.txt --out AFR_vs_notAFR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AMR_indivfam.txt --weir-fst-pop notAMR_indivfam.txt --out AMR_vs_notAMR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop SAS_indivfam.txt --weir-fst-pop notSAS_indivfam.txt --out SAS_vs_notSAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EAS_indivfam.txt --weir-fst-pop notEAS_indivfam.txt --out EAS_vs_notEAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EUR_indivfam.txt --weir-fst-pop notEUR_indivfam.txt --out EUR_vs_notEUR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop OCE_indivfam.txt --weir-fst-pop notOCE_indivfam.txt --out OCE_vs_notOCE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop MDE_indivfam.txt --weir-fst-pop notMDE_indivfam.txt --out MDE_vs_notMDE_fst

vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop AMR_indivfam.txt --out AFR_vs_AMR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop SAS_indivfam.txt --out AFR_vs_SAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop EAS_indivfam.txt --out AFR_vs_EAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop EUR_indivfam.txt --out AFR_vs_EUR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop OCE_indivfam.txt --out AFR_vs_OCE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AFR_indivfam.txt --weir-fst-pop MDE_indivfam.txt --out AFR_vs_MDE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AMR_indivfam.txt --weir-fst-pop SAS_indivfam.txt --out AMR_vs_SAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AMR_indivfam.txt --weir-fst-pop EAS_indivfam.txt --out AMR_vs_EAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AMR_indivfam.txt --weir-fst-pop EUR_indivfam.txt --out AMR_vs_EUR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AMR_indivfam.txt --weir-fst-pop OCE_indivfam.txt --out AMR_vs_OCE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop AMR_indivfam.txt --weir-fst-pop MDE_indivfam.txt --out AMR_vs_MDE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop SAS_indivfam.txt --weir-fst-pop EAS_indivfam.txt --out SAS_vs_EAS_fst
vcftools --gzvcf $vcf2 --weir-fst-pop SAS_indivfam.txt --weir-fst-pop EUR_indivfam.txt --out SAS_vs_EUR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop SAS_indivfam.txt --weir-fst-pop OCE_indivfam.txt --out SAS_vs_OCE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop SAS_indivfam.txt --weir-fst-pop MDE_indivfam.txt --out SAS_vs_MDE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EAS_indivfam.txt --weir-fst-pop EUR_indivfam.txt --out SAS_vs_EUR_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EAS_indivfam.txt --weir-fst-pop OCE_indivfam.txt --out SAS_vs_OCE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EAS_indivfam.txt --weir-fst-pop MDE_indivfam.txt --out SAS_vs_MDE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EUR_indivfam.txt --weir-fst-pop OCE_indivfam.txt --out EUR_vs_OCE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop EUR_indivfam.txt --weir-fst-pop MDE_indivfam.txt --out EUR_vs_MDE_fst
vcftools --gzvcf $vcf2 --weir-fst-pop OCE_indivfam.txt --weir-fst-pop MDE_indivfam.txt --out OCE_vs_MDE_fst
