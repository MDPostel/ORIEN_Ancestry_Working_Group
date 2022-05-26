#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=231G
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=10
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs_HGDP

#---------------------------------------
# DO THE FOLLOWING FOR EACH CHROMOSOME
#--------
i=1 #1-22

ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged_noalts_bedtoolsslop100bp.bed")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
gatk=$(echo "/home1/mpostel/bin/gatk/gatk")

mem=$(echo "230") #memory
t=$(echo "64") #threads
ulimit -n 4001 #files open
date=$(echo March2022)

module load usc
module load openjdk
module load python/3.9.2
module load bzip2
module load bcftools
module load htslib

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/

HGDP=$(echo HGDP_merged_ALLchr.g.vcf.gz)

bcftools view $HGDP -r chr${i} | #bcftools view -m2 -M2 -v snps | #bcftools view -V indels -e 'MAF<0.01' -Oz -o HGDP_merged_chr${i}_snps_maf1.g.vcf.gz --threads $t
bcftools index -t HGDP_merged_chr${i}_snps_maf1.g.vcf.gz --threads $t

echo "subsetting $i by exome bed and dbSNP VCF; filter for DP>20"
file2=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/HGDP_merged_chr${i}_snps_maf1.g.vcf.gz)
name2=$(echo $file | cut -d/ -f9 | cut -d"." -f1)
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/HGDP
mkdir chr${i}
cd chr${i}
bcftools isec -R $bed -p ./ $file2 $dbSNP -Oz --threads $t
mv 0002.vcf.gz ${name2}_subsetbed_dbSNPisec.vcf.gz
bcftools index -t -f ${name2}_subsetbed_dbSNPisec.vcf.gz --threads $t

bcftools sort ${name2}_subsetbed_dbSNPisec.vcf.gz -Oz -o ${name2}_subsetbed_dbSNPisec_sorted.vcf.gz
bcftools index -t -f ${name2}_subsetbed_dbSNPisec_sorted.vcf.gz --threads $t

#---------------------------------------
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset
echo "concatenating"
for i in {1..22}
    do
echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/HGDP/chr${i}/HGDP_merged_chr${i}_snps_maf1_subsetbed_dbSNPisec_sorted.vcf.gz
    done > HGDP_concat_list.txt

bcftools concat -f HGDP_concat_list.txt -Oz -o HGDP_merged_ALLchr_snps_maf1_March2022.vcf.gz --threads $t

echo "imputing"
module load gcc/11.2.0 jdk/17.0.1
/spack/apps2/linux-centos7-x86_64/gcc-11.2.0/jdk-17.0.1-6rjjefs2zw3fedmia6auq5dtifvprffc/bin/java -jar /home1/mpostel/bin/beagle/beagle.25Mar22.4f6.jar gt=HGDP_merged_ALLchr_snps_maf1_March2022.vcf.gz out=HGDP_merged_ALLchr_snps_maf1_impute_March2022 nthreads=20

#---------------------------------------
file2=$(echo HGDP_merged_ALLchr_snps_maf1_impute_March2022.vcf.gz)
bcftools index -t -f $file2 --threads $t
echo "subsetting $file2"
bcftools isec -R $bed $dbSNP $file2 -p ./ -Oz --threads $t

#---------------------------------------
echo "filtering for DP>20"
mv 0003.vcf.gz HGDP_merged_ALLchr_snps_maf1_impute_subsetbybed_dbSNPisec_March2022.vcf.gz
bcftools index -t -f HGDP_merged_ALLchr_snps_maf1_impute_subsetbybed_dbSNPisec_March2022.vcf.gz --threads $t
$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -V HGDP_merged_ALLchr_snps_maf1_impute_subsetbybed_dbSNPisec_March2022.vcf.gz \
   -O HGDP_merged_ALLchr_snps_maf1_impute_subsetbybed_dbSNPisec_March2022.DP.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A DepthPerSampleHC \
   -A Coverage

file3=$(echo HGDP_merged_ALLchr_snps_maf1_impute_subsetbybed_dbSNPisec_March2022.DP.vcf.gz)
prefix2=$(echo $file3 | cut -d"." -f1)
#DP CAN'T BE ANNOTATED ON HGDP?
#bcftools filter --include 'INFO/DP>20' $file3 -Oz -o ${prefix2}_DPgt20.vcf.gz --threads $t
#bcftools index -t -f ${prefix2}_DPgt20.vcf.gz --threads $t


#------
#FILTER FOR PER-COHORT MISSINGNESS <=10%, MAF >=.01 FOR REFERENCES DATASETS, SNPS ONLY
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/
HGDP2=$(echo HGDP_merged_ALLchr_snps_maf1_impute_subsetbybed_dbSNPisec_March2022.DP.vcf.gz)

bcftools view -e 'F_MISSING>0.20' -V indels $HGDP2 -Oz -o HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20.vcf.gz
bcftools index -t -f HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20.vcf.gz --threads $t

HGDP3=$(echo HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20.vcf.gz)


#------
#VQSR

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs

HGDP_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
GATK_ref=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.fasta)

module load bzip2 bcftools htslib vcftools openjdk samtools jdk

HGDP_ref2=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa)
HGDP_ref2=$(echo $HGDP_ref2)

bcftools sort $HGDP3 -Oz -o HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_sorted.vcf.gz -m ${mem}G -T ./; bcftools index -t -f HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_sorted.vcf.gz --threads $t

HGDP_4a=$(echo HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_sorted.vcf.gz)
HGDP_phase3=$(echo ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz)

dbSNP138=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.dbsnp138.vcf)
dbSNP=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz)
omni=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G_omni2.5.hg38.vcf.gz)
hapmap=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/hapmap_3.3.hg38.vcf.gz)

gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)

echo running VariantRecalibrator

$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -V $HGDP_4a \
   -O HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio
bcftools index -t -f HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_sorted_QD.vcf.gz --threads $t
HGDP_4=$(echo HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_sorted_QD.vcf.gz)

$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $HGDP_ref2 \
   --variant $HGDP_4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an FS -an SOR \
   -mode SNP \
   --output HGDP.recal \
   --tranches-file HGDP.tranches \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.5 -tranche 97.0 -tranche 96.5 -tranche 96.0 -tranche 95.5 -tranche 95.0 -tranche 94.5 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.5 -tranche 92.0 -tranche 91.5 -tranche 91.0 -tranche 90.5 -tranche 90.0 \
   --rscript-file HGDP.plots.R

echo running ApplyVQSR

$gatk ApplyVQSR \
   -R $HGDP_ref2 \
   --variant $HGDP_4 \
   --output HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf \
   --tranches-file HGDP.tranches \
   --truth-sensitivity-filter-level 90.0 \
   --recal-file HGDP.recal \
   -mode SNP

#sed -i 's/AD,Number=R/AD,Number=./g' HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf
#grep "AD,Number=" HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf

bgzip -f HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf
bcftools index -t -f HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf.gz --threads $t

##bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP90.00to90.50" || %FILTER=="VQSRTrancheSNP90.50to91.00" || %FILTER=="VQSRTrancheSNP91.00to91.50" || %FILTER=="VQSRTrancheSNP91.50to92.00" || %FILTER=="VQSRTrancheSNP92.00to92.50" || %FILTER=="VQSRTrancheSNP92.50to93.00" || %FILTER=="VQSRTrancheSNP93.00to93.50" || %FILTER=="VQSRTrancheSNP93.50to94.00" || %FILTER=="VQSRTrancheSNP94.00to94.50" || %FILTER=="VQSRTrancheSNP94.50to95.00" || %FILTER=="VQSRTrancheSNP95.00to95.50" || %FILTER=="VQSRTrancheSNP95.50to96.00" || %FILTER=="VQSRTrancheSNP96.00to96.50" || %FILTER=="VQSRTrancheSNP96.50to97.00" || %FILTER=="VQSRTrancheSNP97.00to97.50" || %FILTER=="VQSRTrancheSNP97.50to98.00" || %FILTER=="VQSRTrancheSNP98.00to98.50" || %FILTER=="VQSRTrancheSNP98.50to99.00"' HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSRPASSplusTopTranches.vcf.gz

bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP99.00to99.50" || %FILTER=="VQSRTrancheSNP99.50to99.90" || %FILTER=="VQSRTrancheSNP99.90to99.95" || %FILTER=="VQSRTrancheSNP99.95to100.00+" || %FILTER=="VQSRTrancheSNP99.95to100.00"' HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSR.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSRPASSplusTopTranches.vcf.gz

bcftools index -t -f HGDP_subsetbed_isecdbSNP_noDPfilter_sorted_snps_maf1_missing20_VQSRPASSplusTopTranches.vcf.gz --threads $t


