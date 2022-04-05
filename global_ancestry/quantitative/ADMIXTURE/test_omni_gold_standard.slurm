#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=test_omni

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/affy_test
module load bzip2 bcftools htslib plink

threads=$(echo 40)
list=$(echo gVCFlist.txt)
ulimit -n 4001
mem=$(echo 239)

#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
#file=$(echo ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz) 
file=$(echo ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz)

ref=$(echo /project/davidwcr_264/Resources/references/Homo_sapiens/gencode/GRCh37.p13.genome.fa)
gatk=$(echo "/project/davidwcr_264/Packages/gatk/4.1.0.0/gatk")
dbSNP=$(echo 00-All_chr.vcf.gz)

bcftools annotate --rename-chr ../subset/nochr_to_chr.txt $file -Oz -o omni_chr.vcf.gz --threads $threads
tabix -f omni_chr.vcf.gz

bcftools view --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 omni_chr.vcf.gz -Oz -o omni_chr_autosomes.vcf.gz
tabix -f omni_chr_autosomes.vcf.gz

bcftools +prune --max LD=0.2 omni_chr_autosomes.vcf.gz -Oz -o omni_chr_autosomes_LDprune.vcf.gz

plink --vcf omni_chr_autosomes_LDprune.vcf.gz --make-bed --out omni

bed=$(echo omni.bed)
K=$(echo 20)
/home1/mpostel/bin/admixture/admixture_linux-1.3.0/admixture ${bed} ${K} -j${threads} --cv=10
