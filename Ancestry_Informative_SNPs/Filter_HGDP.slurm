#!/usr/bin/bash
##SBATCH --partition=dtg
#SBATCH --partition=largemem
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs_HGDP

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HGDP REFERENCE POP
#---------------------------
#SET UP

#ref=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa)
ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged.bed")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
gatk=$(echo "/home1/mpostel/bin/gatk/gatk")

mem=$(echo "199") #memory
t=$(echo "600") #threads
ulimit -n 4001 #files open
date=$(echo August2021)

module load usc
module load openjdk
module load python/3.9.2
module load bzip2
module load bcftools
module load htslib

#---------------------------
#GET THE SAMPLES: https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/gVCFs/

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/
#ls /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/HGDP_ref_files > HGDP_gvcf_list.txt

#---------------------------
#COMBINE GVCFS USING BCFTOOLS (CAN'T USE GATK COMBINEGVCFS FOR THESE)

echo "combining gVCFs"
bcftools merge --file-list HGDP_gvcf_list.txt -Oz -o HGDP_merged.g.vcf.gz --gvcf $ref --threads $t
tabix -f HGDP_merged.g.vcf.gz
   # HELPS TO DO ABOVE STEP IN BATCHES/SUBSETS
   # AND THEN BREAK UP BY CHROMOSOME
for i in {1..22}
    do
bcftools view -R chr${i} HGDP_merged.g.vcf.gz -Oz -o HGDP_merged_chr${i}.g.vcf.gz
tabix -f HGDP_merged_chr${i}.g.vcf.gz

#---------------------------
#EXPAND REF BLOCKS

echo "expanding REF blocks"
file=$(echo HGDP_merged_chr${i}.g.vcf.gz)
name=$(echo $file | cut -d"." -f1)
bcftools convert $file --gvcf2vcf -f $ref -Oz -o ${name}_expandREFblocks.g.vcf.gz
tabix -f ${name}_expandREFblocks.g.vcf.gz

#---------------------------
#SUBSET BY EXOME BED AND DBSNP VCF; FILTER FOR DP>20

echo "subsetting by exome bed and dbSNP VCF; filter for DP>20"
file2=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/${name}_expandREFblocks.g.vcf.gz)
name2=$(echo $file | cut -d"/" -f9 | cut -d"." -f1)
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/HGDP
mkdir chr${i}; cd chr${i}
bcftools isec -R $bed -p ./ $file2 $dbSNP --threads $t
mv 0002.vcf ${name2}_subsetbed_dbSNPisec.vcf
bgzip ${name2}_subsetbed_dbSNPisec.vcf
tabix -f ${name2}_subsetbed_dbSNPisec.vcf.gz

HGDP2=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/HGDP/chr${i}/${name2}_subsetbed_dbSNPisec.vcf.gz")
bcftools filter --include 'INFO/DP>20' $HGDP2 -Oz -o ${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz --threads $t
tabix -f ${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz

#---------------------------
#SORT

echo "sorting"
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/
file3=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/HGDP/chr${i}/${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz")

bcftools sort -m ${mem}G -T ./ $file3 -Oz -o ${name2}_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz #--threads $t
tabix -f ${name2}_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz

    done
#---------------------------
#CONCAT

echo "HGDP_merged_chr1_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr2_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr3_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr4_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr5_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr6_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr7_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr8_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr9_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr10_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr11_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr12_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr13_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr14_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr15_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr16_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr17_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr18_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr19_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr20_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr21_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz
HGDP_merged_chr22_expandREFblocks_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz" > HGDP_concat_list.txt

bcftools concat -f HGDP_concat_list.txt -Oz -o HGDP_final.vcf.gz --threads $t
tabix -f HGDP_final.vcf.gz
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
