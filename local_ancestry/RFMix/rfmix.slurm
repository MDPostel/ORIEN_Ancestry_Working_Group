#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=rfmixchr1

threads=$SLURM_CPUS_PER_TASK
module load bzip2 bcftools htslib openjdk
chr=$(echo 1)

export PATH="/home1/mpostel/bin/rfmix:$PATH"

#plink --bfile ${bed} --recode --out isec

echo Activating conda py-popgen environment to run beagle and phase VCFs
#conda init
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda create -n py-popgen python=3.7.7
#conda activate py-popgen
#conda install -c jaredgk -c bioconda py-popgen
##conda install -c bioconda java-jdk
#conda install -c anaconda openjdk
conda activate py-popgen
#module load bzip2 bcftools htslib openjdk

vcf_phase.py --vcf $query_vcf --phase-algorithm beagle --out ORIEN_TCGA_chrfix_phased --out-format vcf.gz
vcf_phase.py --vcf $ref_vcf --phase-algorithm beagle --out HGDP_1000GP_chrfix_phased --out-format vcf.gz

tabix -f *phased.vcf.gz

out=$(echo rfmix_chr${i})

echo Running rfmix
rfmix -f ${query_vcf_phased} -r ${ref_vcf2_phased} -m ${sample_map} -g ${genetic_map} -o ${out} --chromosome=${chr} --n-threads=${threads}
