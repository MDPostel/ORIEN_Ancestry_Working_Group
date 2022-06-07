#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=8G
#SBATCH --cpus-per-task=64
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=rfmixchr1

threads=$SLURM_CPUS_PER_TASK

chr=$(echo 1)


#plink --bfile ${bed} --recode --out isec

echo Activating conda py-popgen environment to run beagle and phase VCFs
#conda create -n py-popgen
conda activate py-popgen
#conda install -c bioconda py-popgen

vcf_phase.py --vcf $query_vcf --phase-algorithm beagle --out ORIEN_TCGA_chrfix_phased --out-format vcf.gz
vcf_phase.py --vcf $ref_vcf --phase-algorithm beagle --out HGDP_1000GP_chrfix_phased --out-format vcf.gz

tabix -f *phased.vcf.gz

out=$(echo rfmix_chr${i})

echo Running rfmix
rfmix -f ${query_vcf_phased} -r ${ref_vcf2_phased} -m ${sample_map} -g ${genetic_map} -o ${out} --chromosome=${chr} --n-threads=${threads} --crf-weight=3
