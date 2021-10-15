#!/usr/bin/bash
#SBATCH --partition=largemem
##SBATCH --partition=oneweek
#SBATCH --mem=250G
#SBATCH --time=07-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
ulimit -n 4001
vcftools=$(echo "/home1/mpostel/bin/vcftools/src/cpp/vcftools")
plink=$(echo "/home1/mpostel/bin/plink-1.90/plink")
t=$(echo 4000) #threads
mem=$(echo 249) #memory

module load bzip2 bcftools htslib vcftools openjdk samtools 

HGDP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/HGDP_final_Sep2021.vcf.gz")
GP1000=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20_re-sorted.vcf.gz")
TCGA=$(echo "/project/davidwcr_666/TCGA/WES/TCGA_final.vcf.gz")

#------
#************************************************ ONLY THING TO CHANGE FOR DIFFERENT ORIEN SAMPLES: ************************************************

date=$(echo August2021)
ORIEN=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted_reheader.vcf.gz")

#***************************************************************************************************************************************************

#------
#FILTER FOR PER-COHORT MISSINGNESS <=10%, MAF >=.01 FOR REFERENCES DATASETS, SNPS ONLY

bcftools view $HGDP -m2 -M2 -v snps | bcftools view -V indels -e 'MAF<0.01' -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz --threads $t
tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz
HGDP2=$(echo "HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz")

bcftools view $ORIEN -m2 -M2 -v snps | bcftools view -V indels -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz --threads $t
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz
ORIEN2=$(echo "ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz")

bcftools view $GP1000 -m2 -M2 -v snps | bcftools view -V indels -e 'MAF<0.01' -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz --threads $t
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz
GP1000_2=$(echo "GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.vcf.gz")

bcftools view $TCGA -m2 -M2 -v snps | bcftools view -V indels -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz --threads $t
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz
TCGA2=$(echo "TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz")
bcftools view -e 'F_MISSING>0.10' -V indels $TCGA2 -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz

tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
HGDP3=$(echo "HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz")
ORIEN3=$(echo "ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz")
GP1000_3=$(echo "GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz")
TCGA3=$(echo "TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz")

#------
#VQSR

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs

ORIEN_ref=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa)
TCGA_ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")
GP1000_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
HGDP_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
GATK_ref=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.fasta)

cat $ORIEN_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa
cat ${ORIEN_ref}.fai | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa.fai
java -jar /project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar CreateSequenceDictionary \
R=/project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa \
O=/project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.dict

cat $TCGA_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa
cat ${TCGA_ref}.fai | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa.fai
java -jar /project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar CreateSequenceDictionary \
R=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa \
O=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.dict
TCGA_ref2=$(echo /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa)

cat $GP1000_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa
cat ${GP1000_ref}.fai | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa.fai
java -jar /project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar CreateSequenceDictionary \
R=/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa \
O=/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.dict
GP1000_ref2=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa)
HGDP_ref2=$(echo $GP1000_ref2)

bcftools view -h $HGDP3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > HGDP_header.txt
sed -i 's/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla_fix/g' HGDP_header.txt
bcftools view -h $GP1000_3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > 1000GP_header.txt
sed -i 's/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla_fix/g' 1000GP_header.txt
bcftools view -h $TCGA3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > TCGA_header.txt
sed -i 's/GRCh38.primary_assembly.genome/GRCh38.primary_assembly.genome_fix/g' TCGA_header.txt
bcftools view -h $ORIEN3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > ORIEN_header.txt
sed -i 's/hs38DH/hs38DH_fix/g' ORIEN_header.txt

bcftools reheader -h HGDP_header.txt -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_reheader.vcf HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
bcftools reheader -h 1000GP_header.txt -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_reheader.vcf GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
bcftools reheader -h TCGA_header.txt -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_reheader.vcf TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gzbcftools reheader -h ORIEN_header.txt -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_reheader.vcf ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz

bcftools sort HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz
bcftools sort GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz
bcftools sort ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz
bcftools sort TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz

HGDP4a=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz)
GP1000_4a=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz)
ORIEN4=$(echo ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz)
TCGA4a=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz)

GP1000_phase3=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf)
GP1000_phase1=$(echo 1000G_phase1.snps.high_confidence.hg38.vcf.gz)
GP1000_phase3=$(echo ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz)
dbSNP138=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.dbsnp138.vcf)
dbSNP=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz)
omni=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G_omni2.5.hg38.vcf.gz)
hapmap=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/hapmap_3.3.hg38.vcf.gz)

gatk=$(echo /home1/mpostel/bin/gatk/gatk)
#gatk=$(echo /project/davidwcr_264/Packages/gatk/4.0.10.1/gatk)
#gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.0.0/gatk)
#gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)
#gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.5.0/gatk)
#gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.4.1/gatk)

#---
#REMOVED PARAMETERS FROM VARIANTRECALIBRATOR FOR ALL 4 DATASETS:
#   --disable-sequence-dictionary-validation \
#   --resource:1000G,known=false,training=true,truth=false,prior=10.0 $GP1000_phase3 \
#   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP \

#---
$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $ORIEN_ref2 \
   --variant $ORIEN4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   --output ORIEN.recal \
   --tranches-file ORIEN.tranches \
   --rscript-file ORIEN.plots.R

#---
$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $TCGA_ref2 \
   -V $TCGA4a \
   -O TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \   
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio

tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted_QD.vcf.gz
TCGA4=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted_QD.vcf.gz)

#---
#REMOVED PARAMETERS:
   #-an QD
   #-an MQRankSum
   #-an ReadPosRankSum
$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $TCGA_ref2 \
   --variant $TCGA4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an MQ -an FS -an SOR \
   -mode SNP \
   --output TCGA.recal \
   --tranches-file TCGA.tranches \
   --rscript-file TCGA.plots.R

#---
#REMOVED PARAMETERS:
   #   -R $HGDP_ref2 \
$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -V $HGDP4a \
   -O HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio

tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted_QD.vcf.gz
HGDP4=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted_QD.vcf.gz)

#---
#REMOVED PARAMETERS:
   #-an MQ
   #-an MQRankSum
   #-an ReadPosRankSum
$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $HGDP_ref2 \
   --variant $HGDP4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an QD -an FS -an SOR \
   -mode SNP \
   --output HGDP.recal \
   --tranches-file HGDP.tranches \
   --rscript-file HGDP.plots.R

#---
#REMOVED PARAMETERS:
#   -R $GP1000_ref2 \
$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -V $GP1000_4a \
   -O GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted_QD.vcf.gz
GP1000_4=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted_QD.vcf.gz)

#---
#REMOVED PARAMETERS:#REMOVED PARAMETERS:
   #-an QD
   #-an MQ
   #-an MQRankSum
   #-an ReadPosRankSum
$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $GP1000_ref2 \
   --variant $GP1000_4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an FS -an SOR \
   -mode SNP \
   --output 1000GP.recal \
   --tranches-file 1000GP.tranches \
   --rscript-file 1000GP.plots.R

#---
$gatk ApplyVQSR \
   -R $ORIEN_ref2 \
   --variant $ORIEN4 \
   --output ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file ORIEN.tranches \
   --recal-file ORIEN.recal \
   -mode SNP

$gatk ApplyVQSR \
   -R $TCGA_ref2 \
   --variant $TCGA4 \
   --output TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file TCGA.tranches \
   --recal-file TCGA.recal \
   -mode SNP

$gatk ApplyVQSR \
   -R $HGDP_ref2 \
   --variant $HGDP4 \
   --output HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file HGDP.tranches \
   --recal-file HGDP.recal \
   -mode SNP

$gatk ApplyVQSR \
   -R $GP1000_ref2 \
   --variant $GP1000_4 \
   --output GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file 1000GP.tranches \
   --recal-file 1000GP.recal \
   -mode SNP

sed -i 's/FORMAT=<ID=AD,Number=R/FORMAT=<ID=AD,Number=./g' *VQSR.vcf
bgzip ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf
bgzip TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf
bgzip HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf
bgzip GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz
tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz

#------
#MERGE:
HGDP5=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz)
ORIEN5=$(echo ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz)
GP1000_5=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz)
TCGA5=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz)

bcftools merge $HGDP5 $GP1000_5 $ORIEN5 $TCGA5 --gvcf $ORIEN_ref2 -Oz -o HGDP_1000GP_ORIEN_TCGA_snps_missing10_merged.vcf.gz --threads $t
tabix -f HGDP_1000GP_ORIEN_TCGA_snps_missing10_merged.vcf.gz

merge=$(echo "HGDP_1000GP_ORIEN_TCGA_snps_missing10_merged.vcf.gz")

#-----------------
#VCFTOOLS VERSION

vcftools --gzvcf $merge --max-missing 0.90 --recode --out merge
# --remove-indels --chr 1-22

bgzip merge.recode.vcf; tabix -f merge.recode.vcf.gz
merge4=$(echo merge.recode.vcf.gz)

bcftools view $merge4 --force-samples -S HGDP_sample_names.txt --threads $t > HGDP_postmerge.vcf
bcftools view $merge4 --force-samples -S ORIEN_sample_names.txt --threads $t > ORIEN_postmerge.vcf
bcftools view $merge4 --force-samples -S 1000GP_sample_names.txt --threads $t > 1000GP_postmerge.vcf
bcftools view $merge4 --force-samples -S TCGA_sample_names.txt --threads $t > TCGA_postmerge.vcf
bgzip HGDP_postmerge.vcf; tabix -f HGDP_postmerge.vcf.gz
bgzip ORIEN_postmerge.vcf; tabix -f ORIEN_postmerge.vcf.gz
bgzip 1000GP_postmerge.vcf; tabix -f 1000GP_postmerge.vcf.gz
bgzip TCGA_postmerge.vcf; tabix -f TCGA_postmerge.vcf.gz

cat HGDP_sample_names.txt > HGDP_1000GP_samples.txt
cat 1000GP_sample_names.txt >> HGDP_1000GP_samples.txt
bcftools view $merge4 --force-samples -S HGDP_1000GP_samples.txt --threads $t -Oz -o HGDP_1000GP.vcf.gz
tabix -f HGDP_1000GP.vcf.gz

bcftools view HGDP_1000GP.vcf.gz --force-samples -S AFR.txt -Oz -o AFR.vcf.gz --threads $t
bcftools view HGDP_1000GP.vcf.gz --force-samples -S EUR.txt -Oz -o EUR.vcf.gz --threads $t
bcftools view HGDP_1000GP.vcf.gz --force-samples -S AMR.txt -Oz -o AMR.vcf.gz --threads $t
bcftools view HGDP_1000GP.vcf.gz --force-samples -S SAS.txt -Oz -o SAS.vcf.gz --threads $t
bcftools view HGDP_1000GP.vcf.gz --force-samples -S EAS.txt -Oz -o EAS.vcf.gz --threads $t
bcftools view HGDP_1000GP.vcf.gz --force-samples -S OCE.txt -Oz -o OCE.vcf.gz --threads $t
bcftools view HGDP_1000GP.vcf.gz --force-samples -S MDE.txt -Oz -o MDE.vcf.gz --threads $t

tabix -f AFR.vcf.gz; tabix -f EUR.vcf.gz; tabix -f AMR.vcf.gz; tabix -f SAS.vcf.gz; tabix -f EAS.vcf.gz; tabix -f OCE.vcf.gz; tabix -f MDE.vcf.gz

vcftools --gzvcf $merge4 --hardy --out pooled_postmerge
vcftools --gzvcf $HGDP5 --hardy --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --hardy --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --hardy --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --hardy --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --hardy --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --hardy --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --hardy --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --hardy --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --hardy --out AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --hardy --out EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --hardy --out AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --hardy --out SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --hardy --out EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --hardy --out OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --hardy --out MDE_postmerge


vcftools --gzvcf $merge4 --missing-indv --out pooled_postmerge
vcftools --gzvcf $HGDP5 --missing-indv --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --missing-indv --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --missing-indv --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --missing-indv --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --missing-indv --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --missing-indv --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --missing-indv --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --missing-indv --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --missing-indv --out AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --missing-indv --out EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --missing-indv --out AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --missing-indv --out SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --missing-indv --out EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --missing-indv --out OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --missing-indv --out MDE_postmerge

vcftools --gzvcf $merge4 --missing-site --out pooled_postmerge
vcftools --gzvcf $HGDP5 --missing-site --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --missing-site --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --missing-site --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --missing-site --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --missing-site --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --missing-site --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --missing-site --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --missing-site --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --missing-site --out AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --missing-site --out EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --missing-site --out AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --missing-site --out SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --missing-site --out EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --missing-site --out OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --missing-site --out MDE_postmerge

vcftools --gzvcf $merge4 --freq --out pooled_postmerge
vcftools --gzvcf $HGDP5 --freq --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --freq --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --freq --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --freq --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --freq --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --freq --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --freq --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --freq --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --freq --out AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --freq --out EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --freq --out AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --freq --out SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --freq --out EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --freq --out OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --freq --out MDE_postmerge

vcftools --gzvcf $merge4 --depth --out pooled_postmerge
vcftools --gzvcf $HGDP5 --depth --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --depth --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --depth --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --depth --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --depth --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --depth --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --depth --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --depth --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --depth --out AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --depth --out EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --depth --out AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --depth --out SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --depth --out EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --depth --out OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --depth --out MDE_postmerge

vcftools --gzvcf $merge4 --site-depth --out pooled_postmerge
vcftools --gzvcf $HGDP5 --site-depth --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --site-depth --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --site-depth --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --site-depth --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --site-depth --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --site-depth --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --site-depth --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --site-depth --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --depth --site-depth AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --depth --site-depth EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --depth --site-depth AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --depth --site-depth SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --depth --site-depth EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --depth --site-depth OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --depth --site-depth MDE_postmerge

vcftools --gzvcf $merge4 --hap-r2 --out pooled_postmerge
vcftools --gzvcf $HGDP5 --hap-r2 --out HGDP_PREmerge
vcftools --gzvcf $GP1000_5 --hap-r2 --out 1000GP_PREmerge
vcftools --gzvcf $TCGA5 --hap-r2 --out TCGA_PREmerge
vcftools --gzvcf $ORIEN5 --hap-r2 --out ORIEN_PREmerge
vcftools --gzvcf HGDP_postmerge.vcf.gz --hap-r2 --out HGDP_postmerge
vcftools --gzvcf ORIEN_postmerge.vcf.gz --hap-r2 --out ORIEN_postmerge
vcftools --gzvcf 1000GP_postmerge.vcf.gz --hap-r2 --out 1000GP_postmerge
vcftools --gzvcf TCGA_postmerge.vcf.gz --hap-r2 --out TCGA_postmerge
vcftools --gzvcf AFR.vcf.gz --hap-r2 --out AFR_postmerge
vcftools --gzvcf EUR.vcf.gz --hap-r2 --out EUR_postmerge
vcftools --gzvcf AMR.vcf.gz --hap-r2 --out AMR_postmerge
vcftools --gzvcf SAS.vcf.gz --hap-r2 --out SAS_postmerge
vcftools --gzvcf EAS.vcf.gz --hap-r2 --out EAS_postmerge
vcftools --gzvcf OCE.vcf.gz --hap-r2 --out OCE_postmerge
vcftools --gzvcf MDE.vcf.gz --hap-r2 --out MDE_postmerge
                                                                                                                                                                                        
