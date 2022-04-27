#!/usr/bin/bash
##SBATCH --partition=largemem
#SBATCH --partition=epyc-64
#SBATCH --mem=230G
#SBATCH --cpus-per-task=64
#SBATCH --time=2-00:00:00
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
bed=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged_noalts_bedtoolsslop100bp.bed)
gatk=$(echo "/project/davidwcr_264/Packages/gatk/4.1.0.0/gatk")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
mem=$(echo 230) #memory
t=$(echo 64) #threads
ulimit -n 4001 #files open

module load usc bzip2 bcftools htslib r 

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
#mkdir 1000GP
cd 1000GP

for i in {1..22}
    do
cd ${dir}/isec/1000GP
#mkdir chr${i}
cd chr${i}

file2=$(echo "${dir}/1000GP.chr${i}_sorted.vcf.gz")
echo "subsetting $file2"
bcftools isec -R $bed $dbSNP $file2 -p ${dir}/isec/1000GP/chr${i} -Oz --threads $t

#---------------------------
#FILTER FOR DP>20

echo "filtering for DP>20 -- chr$i"
mv ${dir}/isec/1000GP/chr${i}/0003.vcf.gz ${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf.gz
bcftools index -t -f ${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf.gz --threads $t

file3=$(echo "${dir}/isec/1000GP/chr${i}/1000GP.chr${i}_sorted_subsetbybed_dbSNPisec.vcf.gz")

cd $dir
echo "filtering $file3 for DP>20"
prefix2=$(echo $file3 | cut -d"/" -f12 | cut -d"." -f1,2)
bcftools filter --include 'DP>20' $file3 -Oz -o ${prefix2}_DPgt20.vcf.gz --threads $t
tabix -f ${prefix2}_DPgt20.vcf.gz
    done
    
#---------------------------
#MERGE AND RE-SORT

echo "merging & re-sorting"
cd $dir
for i in {1..22}; do echo "1000GP.chr${i}_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz"; done > 1000GP_chr_list_to_concat.txt

bcftools concat -f 1000GP_chr_list_to_concat.txt -Oz -o 1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz --threads $t
bcftools index -t -f 1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz --threads $t

file4=$(echo "1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20.vcf.gz")
prefix3=$(echo $file4 | cut -d"." -f1,2)
bcftools sort -m ${mem}G -T $dir $file4 -Oz -o ${prefix3}_re-sorted_March2022.vcf.gz #--threads $t
tabix -f ${prefix3}_re-sorted_March2022.vcf.gz


#------
#IMPUTE, FILTER FOR PER-COHORT MISSINGNESS <=20%, MAF >=.01 FOR REFERENCES DATASETS, SNPS ONLY
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/

GP1000=$(echo ${prefix3}_re-sorted_March2022.vcf.gz)
bcftools view $GP1000 -m2 -M2 -v snps | bcftools view -V indels -e 'MAF<0.01' -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz --threads $t
bcftools index -t -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz --threads $t

echo "imputing"
module load gcc/11.2.0 jdk/17.0.1
/spack/apps2/linux-centos7-x86_64/gcc-11.2.0/jdk-17.0.1-6rjjefs2zw3fedmia6auq5dtifvprffc/bin/java -jar /home1/mpostel/bin/beagle/beagle.25Mar22.4f6.jar gt=GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz out=GP1000_merged_ALLchr_snps_maf1_impute_March2022 nthreads=50

GP1000_2=$(echo "GP1000_merged_ALLchr_snps_maf1_impute_March2022.vcf.gz")
bcftools index -t -f $GP1000_2 --threads $t

bcftools view -e 'F_MISSING>0.20' -V indels $GP1000_2 -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz
bcftools index -t -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz --threads $t

GP1000_3=$(echo "GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz")

#------
#VQSR

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs

GP1000_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")

GATK_ref=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.fasta)

module load bzip2 bcftools htslib vcftools openjdk samtools jdk

GP1000_ref2=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa)
HGDP_ref2=$(echo $GP1000_ref2)

bcftools sort GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz -m ${mem}G -T ./; bcftools index -t -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz --threads $t

GP1000_4a=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz)

GP1000_phase3=$(echo ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz)

dbSNP138=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.dbsnp138.vcf)
dbSNP=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz)

omni=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G_omni2.5.hg38.vcf.gz)
hapmap=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/hapmap_3.3.hg38.vcf.gz)

gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)

echo running VariantRecalibrator

$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -V $GP1000_4a \
   -O GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio
bcftools index -t -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz --threads $t
GP1000_4=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz)

$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $GP1000_ref2 \
   --variant $GP1000_4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an FS -an SOR \
   -mode SNP \
   --output 1000GP.recal \
   --tranches-file 1000GP.tranches \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.5 -tranche 97.0 -tranche 96.5 -tranche 96.0 -tranche 95.5 -tranche 95.0 -tranche 94.5 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.5 -tranche 92.0 -tranche 91.5 -tranche 91.0 -tranche 90.5 -tranche 90.0 \
   --rscript-file 1000GP.plots.R

echo running ApplyVQSR

$gatk ApplyVQSR \
   -R $GP1000_ref2 \
   --variant $GP1000_4 \
   --output GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf \
   --tranches-file 1000GP.tranches \
   --truth-sensitivity-filter-level 90.0 \
   --recal-file 1000GP.recal \
   -mode SNP

#sed -i 's/AD,Number=R/AD,Number=./g' GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf
#grep "AD,Number=" GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf

bgzip -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf
bcftools index -t -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz --threads $t

bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP90.00to90.50" || %FILTER=="VQSRTrancheSNP90.50to91.00" || %FILTER=="VQSRTrancheSNP91.00to91.50" || %FILTER=="VQSRTrancheSNP91.50to92.00" || %FILTER=="VQSRTrancheSNP92.00to92.50" || %FILTER=="VQSRTrancheSNP92.50to93.00" || %FILTER=="VQSRTrancheSNP93.00to93.50" || %FILTER=="VQSRTrancheSNP93.50to94.00" || %FILTER=="VQSRTrancheSNP94.00to94.50" || %FILTER=="VQSRTrancheSNP94.50to95.00" || %FILTER=="VQSRTrancheSNP95.00to95.50" || %FILTER=="VQSRTrancheSNP95.50to96.00" || %FILTER=="VQSRTrancheSNP96.00to96.50" || %FILTER=="VQSRTrancheSNP96.50to97.00" || %FILTER=="VQSRTrancheSNP97.00to97.50" || %FILTER=="VQSRTrancheSNP97.50to98.00" || %FILTER=="VQSRTrancheSNP98.00to98.50" || %FILTER=="VQSRTrancheSNP98.50to99.00"' GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASSplusTopTranches.vcf.gz

#bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP99.00to99.50" || %FILTER=="VQSRTrancheSNP99.50to99.90" || %FILTER=="VQSRTrancheSNP99.90to99.95" || %FILTER=="VQSRTrancheSNP99.95to100.00+" || %FILTER=="VQSRTrancheSNP99.95to100.00"' GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASSplusTopTranches.vcf.gz

bcftools index -t -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASSplusTopTranches.vcf.gz --threads $t
