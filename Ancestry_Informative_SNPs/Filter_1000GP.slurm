#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=1001G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=1000GP


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1000GP = DATABASE OF REFERENCE INDIVIDUALS
#---------------------------
#SET UP

ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged.bed")
gatk=$(echo "/project/davidwcr_264/Packages/gatk/4.1.0.0/gatk")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
mem=$(echo 1000) #memory
t=$(echo 600) #threads
ulimit -n 4001 #files open

module load usc 
module load bzip2 
module load bcftools 
module load htslib 

#wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz*
#ls /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/eigenstrat/1000g_ref_concat.vcf.gz > 1000GP_chr.list
#bcftools concat -f 1000GP_chr.list -Oz -o 1000g_ref_concat.vcf.gz; tabix 1000g_ref_concat.vcf.gz

#ref_1000GP=$(echo "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz")
ref_1000GP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/eigenstrat/1000g_ref_concat.vcf.gz")

#---------------------------
#SPLIT UP CHROMOSOMES AND ADD "CHR" TO NOTATION - EASIER IN TERMS OF MEMORY FOR DOWNSTREAM STEPS

dir=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset")
cd $dir
echo "splitting $ref_1000GP by chromosome"
for i in {1..22}
    do
echo "working on chromosome $i"
bcftools view $ref_1000GP -r ${i} -Oz -o 1000GP.chr${i}_TEMP.vcf.gz --threads $t
tabix 1000GP.chr${i}_TEMP.vcf.gz
bcftools annotate --rename-chrs nochr_to_chr.txt 1000GP.chr${i}_TEMP.vcf.gz -Oz -o 1000GP.chr${i}.vcf.gz --threads $t
tabix 1000GP.chr${i}.vcf.gz
rm 1000GP.chr*_TEMP.vcf.gz*
    done

#---------------------------
#SORT

for i in {1..22}
    do
file=$(echo "${dir}/1000GP.chr${i}.vcf.gz")
prefix=$(echo $file | cut -d"/" -f9 | cut -d"." -f1,2)
echo "sorting $file"
bcftools sort -m ${mem}G -T ./ $file -Oz -o ${prefix}_sorted.vcf.gz #--threads $t
tabix ${prefix}_sorted.vcf.gz
    done

#---------------------------
#INTERSECTION W/ EXOME BED AND DBSNP VCF

cd ${dir}/isec
mkdir 1000GP
cd 1000GP

for i in {1..22}
    do
cd ${dir}/isec/1000GP
mkdir chr${i}
cd chr${i}

file2=$(echo "${dir}/1000GP.chr${i}_sorted.vcf.gz")
echo "subsetting $file2"
bcftools isec -R $bed $dbSNP $file2 -p ${dir}/isec/1000GP/chr${i} --threads $t
    done

#---------------------------
#FILTER FOR DP>20

for i in {1..22}
    do
echo "filtering for DP>20 -- chr$i"
mv ${dir}/isec/1000GP/chr${i}/0003.vcf ${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf
bgzip ${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf
tabix ${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf.gz

file3=$(echo "${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf.gz")

cd $dir
echo "filtering $file3 for DP>20"
prefix2=$(echo $file3 | cut -d"/" -f12 | cut -d"." -f1,2)
bcftools filter --include 'DP>20' $file3 -Oz -o ${prefix2}_DPgt20.vcf.gz --threads $t
tabix ${prefix2}_DPgt20.vcf.gz
    done
    
#---------------------------
#MERGE AND RE-SORT

echo "merging & re-sorting"
cd $dir
find ${dir} -type f -iname "1000GP.chr*_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz" > 1000GP_chr_list_to_concat.txt

bcftools concat -f 1000GP_chr_list_to_concat.txt -Oz -o 1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz --threads $t
tabix 1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz

file4=$(echo "1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz")
prefix3=$(echo $file4 | cut -d"." -f1,2)
bcftools sort -m ${mem}G -T $dir $file4 -Oz -o ${prefix3}_re-sorted.vcf.gz #--threads $t
tabix ${prefix3}_re-sorted.vcf.gz
