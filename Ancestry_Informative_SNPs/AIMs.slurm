#!/usr/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --mem=60G
#SBATCH --cpus-per-task=60
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs_80%geno+VQSR

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
ulimit -n 4001
vcftools=$(echo "/home1/mpostel/bin/vcftools/src/cpp/vcftools")
plink=$(echo "/home1/mpostel/bin/plink-1.90/plink")
t=$(echo 60) #threads
mem=$(echo 60) #memory

module load usc; module load bzip2; module load bcftools; module load htslib; module load vcftools; module load r

HGDP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/HGDP_final_Sep2021.vcf.gz")

GP1000=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20_re-sorted.vcf.gz")

#------
#************************************************ ONLY THING TO CHANGE FOR DIFFERENT ORIEN SAMPLES: ************************************************ 

date=$(echo Jan2022)

ORIEN=$(echo /project/davidwcr_788/ORIEN_Intermember_Project/isec/ORIEN_intermember_Jan2022_ADfix.vcf.gz)

TCGA=$(echo "/project/davidwcr_666/TCGA/WES/TCGA_CombineGVCFs_Oct2021_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz")

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


bcftools view -e 'F_MISSING>0.20' -V indels $HGDP2 -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz
bcftools view -e 'F_MISSING>0.20' -V indels $ORIEN2 -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz
bcftools view -e 'F_MISSING>0.20' -V indels $GP1000_2 -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz
bcftools view -e 'F_MISSING>0.20' -V indels $TCGA2 -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz

tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz

HGDP3=$(echo "HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz")
ORIEN3=$(echo "ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz")
GP1000_3=$(echo "GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz")
TCGA3=$(echo "TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz")

#------
#VQSR

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs

ORIEN_ref=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa)
TCGA_ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")
GP1000_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
HGDP_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
GATK_ref=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.fasta)

module load bzip2 bcftools htslib vcftools openjdk samtools jdk

#cat $ORIEN_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa
#cat ${ORIEN_ref}.fai | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa.fai
#java -jar /project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar CreateSequenceDictionary \
#R=/project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa \
#O=/project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.dict
ORIEN_ref2=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hs38DH_fix.fa)

#cat $TCGA_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa
#cat ${TCGA_ref}.fai | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa.fai
#java -jar /project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar CreateSequenceDictionary \
#R=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa \
#O=/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.dict
TCGA_ref2=$(echo /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa)

#cat $GP1000_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa
#dos2unix /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa
#sed -n 1,24p /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai > /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa.fai
#java -jar /project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar CreateSequenceDictionary \
#R=/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa \
#O=/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.dict
GP1000_ref2=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa)
HGDP_ref2=$(echo $GP1000_ref2)

#bcftools view -h $HGDP3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > HGDP_header.txt
#sed -i 's/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla_fix/g' HGDP_header.txt
#bcftools view -h $GP1000_3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > 1000GP_header.txt
#sed -i 's/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla_fix/g' 1000GP_header.txt
#bcftools view -h $TCGA3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > TCGA_header.txt
#sed -i 's/GRCh38.primary_assembly.genome/GRCh38.primary_assembly.genome_fix/g' TCGA_header.txt
#bcftools view -h $ORIEN3 | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > ORIEN_header.txt
#sed -i 's/hs38DH/hs38DH_fix/g' ORIEN_header.txt


bcftools sort HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz -m ${mem}G -T ./; tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz
bcftools sort GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz -m ${mem}G -T ./; tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz
bcftools sort ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted.vcf.gz -m ${mem}G -T ./; tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted.vcf.gz
bcftools sort TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20.vcf.gz -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted.vcf.gz -m ${mem}G -T ./; tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted.vcf.gz

HGDP4a=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz)
GP1000_4a=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted.vcf.gz)
ORIEN4=$(echo ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted.vcf.gz)
TCGA4a=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted.vcf.gz)


GP1000_phase3=$(echo ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz)
dbSNP138=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.dbsnp138.vcf)
dbSNP=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz)
omni=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G_omni2.5.hg38.vcf.gz)
hapmap=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/hapmap_3.3.hg38.vcf.gz)
gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)


echo running VariantRecalibrator
$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $ORIEN_ref2 \
   --variant $ORIEN4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.5 -tranche 97.0 -tranche 96.5 -tranche 96.0 -tranche 95.5 -tranche 95.0 -tranche 94.5 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.5 -tranche 92.0 -tranche 91.5 -tranche 91.0 -tranche 90.5 -tranche 90.0 \
   --output ORIEN.recal \
   --tranches-file ORIEN.tranches \
   --rscript-file ORIEN.plots.R

$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $TCGA_ref2 \
   -V $TCGA4a \
   -O TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio
 
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted_QD.vcf.gz
TCGA4=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_sorted_QD.vcf.gz)

$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $TCGA_ref2 \
   --variant $TCGA4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an QD -an MQRankSum -an ReadPosRankSum -an MQ -an FS -an SOR \
   -mode SNP \
   --output TCGA.recal \
   --max-gaussians 4 \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.5 -tranche 97.0 -tranche 96.5 -tranche 96.0 -tranche 95.5 -tranche 95.0 -tranche 94.5 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.5 -tranche 92.0 -tranche 91.5 -tranche 91.0 -tranche 90.5 -tranche 90.0 \
   --tranches-file TCGA.tranches \
   --rscript-file TCGA.plots.R

$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -V $HGDP4a \
   -O HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz \
   -A QualByDepth \
   -A AS_QualByDepth \
   -A RMSMappingQuality \
   -A MappingQualityRankSumTest \
   -A ReadPosRankSumTest \
   -A FisherStrand \
   -A StrandOddsRatio

tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz
HGDP4=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz)


$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $HGDP_ref2 \
   --variant $HGDP4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an QD -an FS -an SOR \
   -mode SNP \
   --output HGDP.recal \
   --tranches-file HGDP.tranches \
   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.5 -tranche 97.0 -tranche 96.5 -tranche 96.0 -tranche 95.5 -tranche 95.0 -tranche 94.5 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.5 -tranche 92.0 -tranche 91.5 -tranche 91.0 -tranche 90.5 -tranche 90.0 \
   --rscript-file HGDP.plots.R

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
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_sorted_QD.vcf.gz
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
   -R $ORIEN_ref2 \
   --variant $ORIEN4 \
   --output ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf \
   --tranches-file ORIEN.tranches \
   --truth-sensitivity-filter-level 90.0 \
   --recal-file ORIEN.recal \
   -mode SNP

$gatk ApplyVQSR \
   -R $TCGA_ref2 \
   --variant $TCGA4 \
   --output TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf \
   --tranches-file TCGA.tranches \
   --truth-sensitivity-filter-level 90.0 \
   --recal-file TCGA.recal \
   -mode SNP

$gatk ApplyVQSR \
   -R $HGDP_ref2 \
   --variant $HGDP4 \
   --output HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf \
   --tranches-file HGDP.tranches \
   --truth-sensitivity-filter-level 90.0 \
   --recal-file HGDP.recal \
   -mode SNP

$gatk ApplyVQSR \
   -R $GP1000_ref2 \
   --variant $GP1000_4 \
   --output GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf \
   --tranches-file 1000GP.tranches \
   --truth-sensitivity-filter-level 90.0 \
   --recal-file 1000GP.recal \
   -mode SNP

echo ensure that AD definitions match in headers
sed -i 's/AD,Number=R/AD,Number=./g' *missing20_VQSR.vcf
grep "AD,Number=" *missing20_VQSR.vcf

bgzip -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf --threads $t
bgzip -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf --threads $t
bgzip -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf --threads $t
bgzip -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf --threads $t

#bcftools index -t ... --threads $t  ---> might be faster...
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf.gz
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf.gz
tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz


bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP90.00to90.50" || %FILTER=="VQSRTrancheSNP90.50to91.00" || %FILTER=="VQSRTrancheSNP91.00to91.50" || %FILTER=="VQSRTrancheSNP91.50to92.00" || %FILTER=="VQSRTrancheSNP92.00to92.50" || %FILTER=="VQSRTrancheSNP92.50to93.00" || %FILTER=="VQSRTrancheSNP93.00to93.50" || %FILTER=="VQSRTrancheSNP93.50to94.00" || %FILTER=="VQSRTrancheSNP94.00to94.50" || %FILTER=="VQSRTrancheSNP94.50to95.00" || %FILTER=="VQSRTrancheSNP95.00to95.50" || %FILTER=="VQSRTrancheSNP95.50to96.00" || %FILTER=="VQSRTrancheSNP96.00to96.50" || %FILTER=="VQSRTrancheSNP96.50to97.00" || %FILTER=="VQSRTrancheSNP97.00to97.50" || %FILTER=="VQSRTrancheSNP97.50to98.00" || %FILTER=="VQSRTrancheSNP98.00to98.50" || %FILTER=="VQSRTrancheSNP98.50to99.00"' ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf.gz -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSRPASS.vcf.gz

bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP90.00to90.50" || %FILTER=="VQSRTrancheSNP90.50to91.00" || %FILTER=="VQSRTrancheSNP91.00to91.50" || %FILTER=="VQSRTrancheSNP91.50to92.00" || %FILTER=="VQSRTrancheSNP92.00to92.50" || %FILTER=="VQSRTrancheSNP92.50to93.00" || %FILTER=="VQSRTrancheSNP93.00to93.50" || %FILTER=="VQSRTrancheSNP93.50to94.00" || %FILTER=="VQSRTrancheSNP94.00to94.50" || %FILTER=="VQSRTrancheSNP94.50to95.00" || %FILTER=="VQSRTrancheSNP95.00to95.50" || %FILTER=="VQSRTrancheSNP95.50to96.00" || %FILTER=="VQSRTrancheSNP96.00to96.50" || %FILTER=="VQSRTrancheSNP96.50to97.00" || %FILTER=="VQSRTrancheSNP97.00to97.50" || %FILTER=="VQSRTrancheSNP97.50to98.00" || %FILTER=="VQSRTrancheSNP98.00to98.50" || %FILTER=="VQSRTrancheSNP98.50to99.00"' TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSR.vcf.gz -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSRPASS.vcf.gz

bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP90.00to90.50" || %FILTER=="VQSRTrancheSNP90.50to91.00" || %FILTER=="VQSRTrancheSNP91.00to91.50" || %FILTER=="VQSRTrancheSNP91.50to92.00" || %FILTER=="VQSRTrancheSNP92.00to92.50" || %FILTER=="VQSRTrancheSNP92.50to93.00" || %FILTER=="VQSRTrancheSNP93.00to93.50" || %FILTER=="VQSRTrancheSNP93.50to94.00" || %FILTER=="VQSRTrancheSNP94.00to94.50" || %FILTER=="VQSRTrancheSNP94.50to95.00" || %FILTER=="VQSRTrancheSNP95.00to95.50" || %FILTER=="VQSRTrancheSNP95.50to96.00" || %FILTER=="VQSRTrancheSNP96.00to96.50" || %FILTER=="VQSRTrancheSNP96.50to97.00" || %FILTER=="VQSRTrancheSNP97.00to97.50" || %FILTER=="VQSRTrancheSNP97.50to98.00" || %FILTER=="VQSRTrancheSNP98.00to98.50" || %FILTER=="VQSRTrancheSNP98.50to99.00"' HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASS.vcf.gz

bcftools view -i '%FILTER=="PASS" || %FILTER=="VQSRTrancheSNP90.00to90.50" || %FILTER=="VQSRTrancheSNP90.50to91.00" || %FILTER=="VQSRTrancheSNP91.00to91.50" || %FILTER=="VQSRTrancheSNP91.50to92.00" || %FILTER=="VQSRTrancheSNP92.00to92.50" || %FILTER=="VQSRTrancheSNP92.50to93.00" || %FILTER=="VQSRTrancheSNP93.00to93.50" || %FILTER=="VQSRTrancheSNP93.50to94.00" || %FILTER=="VQSRTrancheSNP94.00to94.50" || %FILTER=="VQSRTrancheSNP94.50to95.00" || %FILTER=="VQSRTrancheSNP95.00to95.50" || %FILTER=="VQSRTrancheSNP95.50to96.00" || %FILTER=="VQSRTrancheSNP96.00to96.50" || %FILTER=="VQSRTrancheSNP96.50to97.00" || %FILTER=="VQSRTrancheSNP97.00to97.50" || %FILTER=="VQSRTrancheSNP97.50to98.00" || %FILTER=="VQSRTrancheSNP98.00to98.50" || %FILTER=="VQSRTrancheSNP98.50to99.00"' GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASS.vcf.gz

#could try gatk FilterVariantTranches

tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSRPASS.vcf.gz
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSRPASS.vcf.gz
tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASS.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASS.vcf.gz


echo "CHECK ABOVE FOR ERRORS BEFORE ISEC"

#------
#ISEC -- INTERSECTION OF 4 DATASETS

HGDP5=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASS.vcf.gz)
ORIEN5=$(echo ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSRPASS.vcf.gz)
GP1000_5=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSRPASS.vcf.gz)
TCGA5=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing20_VQSRPASS.vcf.gz)

echo running isec
bcftools isec $HGDP5 $GP1000_5 $ORIEN5 $TCGA5 -n=4 -Oz -o HGDP_1000GP_ORIEN_TCGA_snps_missing20_VQSRPASS_isec.vcf.gz -p /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs

#ISEC DOES NOT END UP NAMING THE OUTPUT CORRECTLY...
echo ensure that AD definitions match in headers
zgrep "AD,Number=" 0000.vcf.gz
zgrep "AD,Number=" 0001.vcf.gz
zgrep "AD,Number=" 0002.vcf.gz
zgrep "AD,Number=" 0003.vcf.gz

bcftools view -h 0000.vcf.gz > 0000.vcf.gz_header.txt
sed -i 's/AD,Number=R/AD,Number=./g' 0000.vcf.gz_header.txt
bcftools reheader -h 0000.vcf.gz_header.txt -o 0000_reheader.vcf.gz 0000.vcf.gz
tabix -f 0000_reheader.vcf.gz
bcftools view -h 0002.vcf.gz > 0002.vcf.gz_header.txt
sed -i 's/AD,Number=R/AD,Number=./g' 0002.vcf.gz_header.txt
bcftools reheader -h 0002.vcf.gz_header.txt -o 0002_reheader.vcf.gz 0002.vcf.gz
tabix -f 0002_reheader.vcf.gz
bcftools view -h 0003.vcf.gz > 0003.vcf.gz_header.txt
sed -i 's/AD,Number=R/AD,Number=./g' 0003.vcf.gz_header.txt
bcftools reheader -h 0003.vcf.gz_header.txt -o 0003_reheader.vcf.gz 0003.vcf.gz
tabix -f 0003_reheader.vcf.gz

bcftools merge 0000_reheader.vcf.gz 0001.vcf.gz 0002_reheader.vcf.gz 0003_reheader.vcf.gz -Oz -o HGDP_1000GP_ORIEN_TCGA_snps_missing20_VQSRPASS_isec.vcf.gz --threads $t 
tabix -f HGDP_1000GP_ORIEN_TCGA_snps_missing20_VQSRPASS_isec.vcf.gz

isec=$(echo HGDP_1000GP_ORIEN_TCGA_snps_missing20_VQSRPASS_isec.vcf.gz)

#TiTv
echo TiTv ORIEN
bcftools stats $ORIEN5 | grep "TSTV"
echo TiTv TCGA
bcftools stats $TCGA5 | grep "TSTV"
echo HGDP
bcftools stats $HGDP5 | grep "TSTV"
echo TiTv 1000GP
bcftools stats $GP1000_5 | grep "TSTV"
echo TiTv isec
bcftools stats $isec | grep "TSTV"

# determine whether there is cryptic relatedness
~/bin/plink-1.90/plink --vcf $isec --double-id --make-bed --out isec
king -b isec.bed --related --degree 2
#------

bcftools view $isec --force-samples -S HGDP_sample_names.txt --threads $t > HGDP_post_missing20_VQSRPASS.vcf
bcftools view $isec --force-samples -S ORIEN_sample_names.txt --threads $t > ORIEN_post_missing20_VQSRPASS.vcf
bcftools view $isec --force-samples -S 1000GP_sample_names.txt --threads $t > 1000GP_post_missing20_VQSRPASS.vcf
bcftools view $isec --force-samples -S TCGA_sample_names.txt --threads $t > TCGA_post_missing20_VQSRPASS.vcf

bgzip -f HGDP_post_missing20_VQSRPASS.vcf; tabix -f HGDP_post_missing20_VQSRPASS.vcf.gz
bgzip -f ORIEN_post_missing20_VQSRPASS.vcf; tabix -f ORIEN_post_missing20_VQSRPASS.vcf.gz
bgzip -f 1000GP_post_missing20_VQSRPASS.vcf; tabix -f 1000GP_post_missing20_VQSRPASS.vcf.gz
bgzip -f TCGA_post_missing20_VQSRPASS.vcf; tabix -f TCGA_post_missing20_VQSRPASS.vcf.gz

cat HGDP_sample_names.txt > HGDP_1000GP_samples.txt
cat 1000GP_sample_names.txt >> HGDP_1000GP_samples.txt
bcftools view $isec --force-samples -S HGDP_1000GP_samples.txt --threads $t -Oz -o HGDP_1000GP_missing20_VQSRPASS.vcf.gz
tabix -f HGDP_1000GP_missing20_VQSRPASS.vcf.gz

#7 samples have to look up on cellosaurus
#echo HGDP00521 >> EUR.txt; echo HGDP00665 >> EUR.txt; echo HGDP00778 >> EAS.txt; echo HGDP01307 >> EAS.txt; echo HGDP00927 >> AFR.txt; echo HGDP01284 >> AFR.txt; echo HGDP00998 >> AMR.txt
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S AFR.txt -Oz -o AFR_missing20_VQSRPASS.vcf.gz --threads $t
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S EUR.txt -Oz -o EUR_missing20_VQSRPASS.vcf.gz --threads $t
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S AMR.txt -Oz -o AMR_missing20_VQSRPASS.vcf.gz --threads $t
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S SAS.txt -Oz -o SAS_missing20_VQSRPASS.vcf.gz --threads $t
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S EAS.txt -Oz -o EAS_missing20_VQSRPASS.vcf.gz --threads $t
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S OCE.txt -Oz -o OCE_missing20_VQSRPASS.vcf.gz --threads $t
bcftools view HGDP_1000GP_missing20_VQSRPASS.vcf.gz --force-samples -S MDE.txt -Oz -o MDE_missing20_VQSRPASS.vcf.gz --threads $t
tabix -f AFR_missing20_VQSRPASS.vcf.gz; tabix -f EUR_missing20_VQSRPASS.vcf.gz; tabix -f AMR_missing20_VQSRPASS.vcf.gz; tabix -f SAS_missing20_VQSRPASS.vcf.gz; tabix -f EAS_missing20_VQSRPASS.vcf.gz; tabix -f OCE_missing20_VQSRPASS.vcf.gz; tabix -f MDE_missing20_VQSRPASS.vcf.gz

#Fst
echo Fst stats
vcftools --gzvcf HGDP_1000GP_missing20_VQSRPASS.vcf.gz --weir-fst-pop AFR.txt --weir-fst-pop AMR.txt --weir-fst-pop SAS.txt --weir-fst-pop EAS.txt --weir-fst-pop EUR.txt --weir-fst-pop OCE.txt --weir-fst-pop MDE.txt --out overall_fsts_missing20_VQSRPASS

vcftools --gzvcf HGDP_1000GP_missing20_VQSRPASS.vcf.gz --weir-fst-pop AFR.txt --weir-fst-pop AMR.txt --weir-fst-pop SAS.txt --weir-fst-pop EAS.txt --weir-fst-pop EUR.txt --out overall_fsts_5trunks_no_OCE_no_MDE_missing20_VQSRPASS

vcftools --gzvcf $isec --hardy --out pooled_post_missing20_VQSRPASS
vcftools --gzvcf $HGDP5 --hardy --out HGDP_pre_missing20_VQSRPASS
vcftools --gzvcf $GP1000_5 --hardy --out 1000GP_pre_missing20_VQSRPASS
vcftools --gzvcf $TCGA5 --hardy --out TCGA_pre_missing20_VQSRPASS
vcftools --gzvcf $ORIEN5 --hardy --out ORIEN_pre_missing20_VQSRPASS
vcftools --gzvcf HGDP_post_missing20_VQSRPASS.vcf.gz --hardy --out HGDP_post_missing20_VQSRPASS
vcftools --gzvcf ORIEN_post_missing20_VQSRPASS.vcf.gz --hardy --out ORIEN_post_missing20_VQSRPASS
vcftools --gzvcf 1000GP_post_missing20_VQSRPASS.vcf.gz --hardy --out 1000GP_post_missing20_VQSRPASS
vcftools --gzvcf TCGA_post_missing20_VQSRPASS.vcf.gz --hardy --out TCGA_post_missing20_VQSRPASS
vcftools --gzvcf AFR_missing20_VQSRPASS.vcf.gz --hardy --out AFR_post_missing20_VQSRPASS
vcftools --gzvcf EUR_missing20_VQSRPASS.vcf.gz --hardy --out EUR_post_missing20_VQSRPASS
vcftools --gzvcf AMR_missing20_VQSRPASS.vcf.gz --hardy --out AMR_post_missing20_VQSRPASS
vcftools --gzvcf SAS_missing20_VQSRPASS.vcf.gz --hardy --out SAS_post_missing20_VQSRPASS
vcftools --gzvcf EAS_missing20_VQSRPASS.vcf.gz --hardy --out EAS_post_missing20_VQSRPASS
vcftools --gzvcf OCE_missing20_VQSRPASS.vcf.gz --hardy --out OCE_post_missing20_VQSRPASS
vcftools --gzvcf MDE_missing20_VQSRPASS.vcf.gz --hardy --out MDE_post_missing20_VQSRPASS

vcftools --gzvcf $isec --missing-site --out pooled_post_missing20_VQSRPASS
vcftools --gzvcf $HGDP5 --missing-site --out HGDP_pre_missing20_VQSRPASS
vcftools --gzvcf $GP1000_5 --missing-site --out 1000GP_pre_missing20_VQSRPASS
vcftools --gzvcf $TCGA5 --missing-site --out TCGA_pre_missing20_VQSRPASS
vcftools --gzvcf $ORIEN5 --missing-site --out ORIEN_pre_missing20_VQSRPASS
vcftools --gzvcf HGDP_post_missing20_VQSRPASS.vcf.gz --missing-site --out HGDP_post_missing20_VQSRPASS
vcftools --gzvcf ORIEN_post_missing20_VQSRPASS.vcf.gz --missing-site --out ORIEN_post_missing20_VQSRPASS
vcftools --gzvcf 1000GP_post_missing20_VQSRPASS.vcf.gz --missing-site --out 1000GP_post_missing20_VQSRPASS
vcftools --gzvcf TCGA_post_missing20_VQSRPASS.vcf.gz --missing-site --out TCGA_post_missing20_VQSRPASS
vcftools --gzvcf AFR_missing20_VQSRPASS.vcf.gz --missing-site --out AFR_post_missing20_VQSRPASS
vcftools --gzvcf EUR_missing20_VQSRPASS.vcf.gz --missing-site --out EUR_post_missing20_VQSRPASS
vcftools --gzvcf AMR_missing20_VQSRPASS.vcf.gz --missing-site --out AMR_post_missing20_VQSRPASS
vcftools --gzvcf SAS_missing20_VQSRPASS.vcf.gz --missing-site --out SAS_post_missing20_VQSRPASS
vcftools --gzvcf EAS_missing20_VQSRPASS.vcf.gz --missing-site --out EAS_post_missing20_VQSRPASS
vcftools --gzvcf OCE_missing20_VQSRPASS.vcf.gz --missing-site --out OCE_post_missing20_VQSRPASS
vcftools --gzvcf MDE_missing20_VQSRPASS.vcf.gz --missing-site --out MDE_post_missing20_VQSRPASS


vcftools --gzvcf $isec --freq --out pooled_post_missing20_VQSRPASS
vcftools --gzvcf $HGDP5 --freq --out HGDP_pre_missing20_VQSRPASS
vcftools --gzvcf $GP1000_5 --freq --out 1000GP_pre_missing20_VQSRPASS
vcftools --gzvcf $TCGA5 --freq --out TCGA_pre_missing20_VQSRPASS
vcftools --gzvcf $ORIEN5 --freq --out ORIEN_pre_missing20_VQSRPASS
vcftools --gzvcf HGDP_post_missing20_VQSRPASS.vcf.gz --freq --out HGDP_post_missing20_VQSRPASS
vcftools --gzvcf ORIEN_post_missing20_VQSRPASS.vcf.gz --freq --out ORIEN_post_missing20_VQSRPASS
vcftools --gzvcf 1000GP_post_missing20_VQSRPASS.vcf.gz --freq --out 1000GP_post_missing20_VQSRPASS
vcftools --gzvcf TCGA_post_missing20_VQSRPASS.vcf.gz --freq --out TCGA_post_missing20_VQSRPASS
vcftools --gzvcf AFR_missing20_VQSRPASS.vcf.gz --freq --out AFR_post_missing20_VQSRPASS
vcftools --gzvcf EUR_missing20_VQSRPASS.vcf.gz --freq --out EUR_post_missing20_VQSRPASS
vcftools --gzvcf AMR_missing20_VQSRPASS.vcf.gz --freq --out AMR_post_missing20_VQSRPASS
vcftools --gzvcf SAS_missing20_VQSRPASS.vcf.gz --freq --out SAS_post_missing20_VQSRPASS
vcftools --gzvcf EAS_missing20_VQSRPASS.vcf.gz --freq --out EAS_post_missing20_VQSRPASS
vcftools --gzvcf OCE_missing20_VQSRPASS.vcf.gz --freq --out OCE_post_missing20_VQSRPASS
vcftools --gzvcf MDE_missing20_VQSRPASS.vcf.gz --freq --out MDE_post_missing20_VQSRPASS

vcftools --gzvcf $isec --depth --out pooled_post_missing20_VQSRPASS
vcftools --gzvcf $HGDP5 --depth --out HGDP_pre_missing20_VQSRPASS
vcftools --gzvcf $GP1000_5 --depth --out 1000GP_pre_missing20_VQSRPASS
vcftools --gzvcf $TCGA5 --depth --out TCGA_pre_missing20_VQSRPASS
vcftools --gzvcf $ORIEN5 --depth --out ORIEN_pre_missing20_VQSRPASS
vcftools --gzvcf HGDP_post_missing20_VQSRPASS.vcf.gz --depth --out HGDP_post_missing20_VQSRPASS
vcftools --gzvcf ORIEN_post_missing20_VQSRPASS.vcf.gz --depth --out ORIEN_post_missing20_VQSRPASS
vcftools --gzvcf 1000GP_post_missing20_VQSRPASS.vcf.gz --depth --out 1000GP_post_missing20_VQSRPASS
vcftools --gzvcf TCGA_post_missing20_VQSRPASS.vcf.gz --depth --out TCGA_post_missing20_VQSRPASS
vcftools --gzvcf AFR_missing20_VQSRPASS.vcf.gz --depth --out AFR_post_missing20_VQSRPASS
vcftools --gzvcf EUR_missing20_VQSRPASS.vcf.gz --depth --out EUR_post_missing20_VQSRPASS
vcftools --gzvcf AMR_missing20_VQSRPASS.vcf.gz --depth --out AMR_post_missing20_VQSRPASS
vcftools --gzvcf SAS_missing20_VQSRPASS.vcf.gz --depth --out SAS_post_missing20_VQSRPASS
vcftools --gzvcf EAS_missing20_VQSRPASS.vcf.gz --depth --out EAS_post_missing20_VQSRPASS
vcftools --gzvcf OCE_missing20_VQSRPASS.vcf.gz --depth --out OCE_post_missing20_VQSRPASS
vcftools --gzvcf MDE_missing20_VQSRPASS.vcf.gz --depth --out MDE_post_missing20_VQSRPASS

vcftools --gzvcf $isec --hap-r2 --out pooled_post_missing20_VQSRPASS
vcftools --gzvcf $HGDP5 --hap-r2 --out HGDP_pre_missing20_VQSRPASS
vcftools --gzvcf $GP1000_5 --hap-r2 --out 1000GP_pre_missing20_VQSRPASS
vcftools --gzvcf $TCGA5 --hap-r2 --out TCGA_pre_missing20_VQSRPASS
vcftools --gzvcf $ORIEN5 --hap-r2 --out ORIEN_pre_missing20_VQSRPASS
vcftools --gzvcf HGDP_post_missing20_VQSRPASS.vcf.gz --hap-r2 --out HGDP_post_missing20_VQSRPASS
vcftools --gzvcf ORIEN_post_missing20_VQSRPASS.vcf.gz --hap-r2 --out ORIEN_post_missing20_VQSRPASS
vcftools --gzvcf 1000GP_post_missing20_VQSRPASS.vcf.gz --hap-r2 --out 1000GP_post_missing20_VQSRPASS
vcftools --gzvcf TCGA_post_missing20_VQSRPASS.vcf.gz --hap-r2 --out TCGA_post_missing20_VQSRPASS
vcftools --gzvcf AFR_missing20_VQSRPASS.vcf.gz --hap-r2 --out AFR_post_missing20_VQSRPASS
vcftools --gzvcf EUR_missing20_VQSRPASS.vcf.gz --hap-r2 --out EUR_post_missing20_VQSRPASS
vcftools --gzvcf AMR_missing20_VQSRPASS.vcf.gz --hap-r2 --out AMR_post_missing20_VQSRPASS
vcftools --gzvcf SAS_missing20_VQSRPASS.vcf.gz --hap-r2 --out SAS_post_missing20_VQSRPASS
vcftools --gzvcf EAS_missing20_VQSRPASS.vcf.gz --hap-r2 --out EAS_post_missing20_VQSRPASS
vcftools --gzvcf OCE_missing20_VQSRPASS.vcf.gz --hap-r2 --out OCE_post_missing20_VQSRPASS
vcftools --gzvcf MDE_missing20_VQSRPASS.vcf.gz --hap-r2 --out MDE_post_missing20_VQSRPASS
