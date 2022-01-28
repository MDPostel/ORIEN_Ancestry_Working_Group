#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=500G
#SBATCH --time=07-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
ulimit -n 4001
plink=$(echo "/home1/mpostel/bin/plink-1.90/plink")
t=$(echo 4000) #threads
mem=$(echo 499) #memory

module load usc; module load bzip2; module load bcftools; module load htslib; module load vcftools; module load r

#IF YOU GET AN ERROR REGARDING AD DEFINITIONS NOT MATCHING AND THIS IMPEDING MERGE -- CHANGE THE AD DEFINITION (##FORMAT=<ID=AD...) IN THE HEADER FROM: "Number=R" TO "Number=."

HGDP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/HGDP_final_Sep2021.vcf.gz")

GP1000=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/1000GP.ALLchr_sorted_subsetbybed_dbSNPisec_DPgt20_re-sorted.vcf.gz")

#------
#************************************************ ONLY THING TO CHANGE FOR DIFFERENT ORIEN SAMPLES: ************************************************ 

date=$(echo Jan2022)
ORIEN=$(echo "/project/davidwcr_788/ORIEN_Intermember_Project/isec/ORIEN_intermember_${date}.vcf.gz")
TCGA=$(echo "/project/davidwcr_666/TCGA/WES/TCGA_CombineGVCFs_Oct2021_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz")


#***************************************************************************************************************************************************
#------
#FILTER FOR PER-COHORT MISSINGNESS <=10%, MAF >=.01 FOR REFERENCES DATASETS, SNPS ONLY

#gunzip $ORIEN
#sed -i 's/AD,Number=R/AD,Number=./g' /project/davidwcr_788/ORIEN_Intermember_Project/isec/ORIEN_intermember_Jan2022.vcf
#bgzip /project/davidwcr_788/ORIEN_Intermember_Project/isec/ORIEN_intermember_Jan2022.vcf
#tabix -f $ORIEN

bcftools view $HGDP -m2 -M2 -v snps | bcftools view -V indels -e 'MAF<0.01' -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.recode.vcf.gz --threads $t
tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.recode.vcf.gz
HGDP2=$(echo "HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.recode.vcf.gz")

bcftools view $ORIEN -m2 -M2 -v snps | bcftools view -V indels -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz --threads $t
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz
ORIEN2=$(echo "ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz")

bcftools view $GP1000 -m2 -M2 -v snps | bcftools view -V indels -e 'MAF<0.01' -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.recode.vcf.gz --threads $t
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.recode.vcf.gz
GP1000_2=$(echo "GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1.recode.vcf.gz")

bcftools view $TCGA -m2 -M2 -v snps | bcftools view -V indels -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz --threads $t
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz
TCGA2=$(echo "TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps.vcf.gz")

bcftools view -e 'F_MISSING>0.10' -V indels $HGDP2 -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
bcftools view -e 'F_MISSING>0.10' -V indels $ORIEN2 -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz
bcftools view -e 'F_MISSING>0.10' -V indels $GP1000_2 -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
bcftools view -e 'F_MISSING>0.10' -V indels $TCGA2 -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz

tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz

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

#bcftools reheader -h HGDP_header.txt -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_reheader.vcf HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
#bcftools reheader -h 1000GP_header.txt -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_reheader.vcf GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz
#bcftools reheader -h TCGA_header.txt -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_reheader.vcf TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz
#bcftools reheader -h ORIEN_header.txt -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_reheader.vcf ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz

#bcftools sort HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz
#bcftools sort GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz
#bcftools sort ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz
#bcftools sort TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10.vcf.gz -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz -m ${mem}G -T ./; tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz

HGDP4a=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz)
GP1000_4a=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_sorted.vcf.gz)
ORIEN4=$(echo ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz)
TCGA4a=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz)


GP1000_phase3=$(echo ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz)
dbSNP138=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.dbsnp138.vcf)
dbSNP=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz)
omni=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G_omni2.5.hg38.vcf.gz)
hapmap=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/hapmap_3.3.hg38.vcf.gz)
gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)

#gunzip $ORIEN4
#sed -i 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf
#bgzip ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf
#tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_sorted.vcf.gz


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

$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
   -R $TCGA_ref2 \
   --variant $TCGA4 \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
   -an QD -an MQRankSum -an ReadPosRankSum -an MQ -an FS -an SOR \
   -mode SNP \
   --output TCGA.recal \
   --max-gaussians 4 \
   --tranches-file TCGA.tranches \
   --rscript-file TCGA.plots.R

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

sed -i 's/AD,Number=R/AD,Number=./g' ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf

rm ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz
bgzip ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf
bgzip TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf
bgzip HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf
bgzip GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf

tabix -f ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz
tabix -f TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz
tabix -f HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz
tabix -f GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz

bcftools view -i '%FILTER=="PASS"' HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz -Oz -o HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR_PASS.vcf.gz
bcftools view -i '%FILTER=="PASS"' ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz -Oz -o ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR_PASS.vcf.gz
bcftools view -i '%FILTER=="PASS"' GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR.vcf.gz -Oz -o GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR_PASS.vcf.gz
bcftools view -i '%FILTER=="PASS"' TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR.vcf.gz -Oz -o TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR_PASS.vcf.gz
tabix -f *PASS.vcf.gz


echo "CHECK ABOVE FOR ERRORS BEFORE ISEC"

#------
#ISEC -- INTERSECTION OF 4 DATASETS

HGDP5=$(echo HGDP_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR_PASS.vcf.gz)
ORIEN5=$(echo ORIEN_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR_PASS.vcf.gz)
GP1000_5=$(echo GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing10_VQSR_PASS.vcf.gz)
TCGA5=$(echo TCGA_subsetbed_isecdbSNP_DPgt20_sorted_snps_missing10_VQSR_PASS.vcf.gz)


bcftools isec $HGDP5 $GP1000_5 $ORIEN5 $TCGA5 -n=4 -Oz -o HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz -p /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs
tabix -f HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz

#AD definitions don't match; reheader
bcftools view -h 0003.vcf.gz -o newheader
sed -i 's/AD,Number=R/AD,Number=./g' newheader
bcftools reheader -h newheader -o 0003_reheader.vcf.gz 0003.vcf.gz
tabix -f 0003_reheader.vcf.gz
rm HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz
bcftools merge 0000.vcf.gz 0001.vcf.gz 0002.vcf.gz 0003_reheader.vcf.gz -Oz -o HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz --threads $t 
tabix -f HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz


isec=$(echo HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz)

#TsTv
bcftools stats $ORIEN5 | grep "TSTV"
bcftools stats $TCGA5 | grep "TSTV"
bcftools stats $HGDP5 | grep "TSTV"
bcftools stats $GP1000_5 | grep "TSTV"
bcftools stats $isec | grep "TSTV"

# determine whether there is cryptic relatedness
~/bin/plink-1.90/plink --vcf $isec --double-id --make-bed --out isec
king -b isec.bed --related --degree 2

#-----------------
#VCFTOOLS VERSION
module load bzip2 bcftools htslib vcftools

rm TCGA.vcf.gz ORIEN.vcf.gz HGDP.vcf.gz GP1000.vcf.gz
mv 0000.vcf.gz HGDP.vcf.gz
mv 0001.vcf.gz GP1000.vcf.gz
mv 0002.vcf.gz ORIEN.vcf.gz
mv 0003_reheader.vcf.gz TCGA.vcf.gz
tabix -f TCGA.vcf.gz ORIEN.vcf.gz HGDP.vcf.gz GP1000.vcf.gz

bcftools merge HGDP.vcf.gz GP1000.vcf.gz --threads $t -Oz -o HGDP_1000GP.vcf.gz
tabix -f HGDP_1000GP.vcf.gz

echo "..."
echo now making vcfs for continental ancestries
echo making AFR
bcftools view HGDP_1000GP.vcf.gz --force-samples -S AFR.txt -Oz -o AFR.vcf.gz --threads $t
echo making EUR
bcftools view HGDP_1000GP.vcf.gz --force-samples -S EUR.txt -Oz -o EUR.vcf.gz --threads $t
echo making AMR
bcftools view HGDP_1000GP.vcf.gz --force-samples -S AMR.txt -Oz -o AMR.vcf.gz --threads $t
echo making SAS
bcftools view HGDP_1000GP.vcf.gz --force-samples -S SAS.txt -Oz -o SAS.vcf.gz --threads $t
echo making EAS
bcftools view HGDP_1000GP.vcf.gz --force-samples -S EAS.txt -Oz -o EAS.vcf.gz --threads $t
echo making OCE
bcftools view HGDP_1000GP.vcf.gz --force-samples -S OCE.txt -Oz -o OCE.vcf.gz --threads $t
echo making MDE
bcftools view HGDP_1000GP.vcf.gz --force-samples -S MDE.txt -Oz -o MDE.vcf.gz --threads $t
tabix -f AFR.vcf.gz; tabix -f EUR.vcf.gz; tabix -f AMR.vcf.gz; tabix -f SAS.vcf.gz; tabix -f EAS.vcf.gz; tabix -f OCE.vcf.gz; tabix -f MDE.vcf.gz

sbatch Fst.slurm
echo "..."
echo gathering HWE stats
vcftools --gzvcf $isec --hardy --out pooled
vcftools --gzvcf HGDP.vcf.gz --hardy --out HGDP
vcftools --gzvcf ORIEN.vcf.gz --hardy --out ORIEN
vcftools --gzvcf 1000GP.vcf.gz --hardy --out 1000GP
vcftools --gzvcf TCGA.vcf.gz --hardy --out TCGA
vcftools --gzvcf AFR.vcf.gz --hardy --out AFR
vcftools --gzvcf EUR.vcf.gz --hardy --out EUR
vcftools --gzvcf AMR.vcf.gz --hardy --out AMR
vcftools --gzvcf SAS.vcf.gz --hardy --out SAS
vcftools --gzvcf EAS.vcf.gz --hardy --out EAS
vcftools --gzvcf OCE.vcf.gz --hardy --out OCE
vcftools --gzvcf MDE.vcf.gz --hardy --out MDE
echo "..."
echo gathering missing site stats
vcftools --gzvcf $isec --missing-site --out pooled
vcftools --gzvcf HGDP.vcf.gz --missing-site --out HGDP
vcftools --gzvcf ORIEN.vcf.gz --missing-site --out ORIEN
vcftools --gzvcf 1000GP.vcf.gz --missing-site --out 1000GP
vcftools --gzvcf TCGA.vcf.gz --missing-site --out TCGA
vcftools --gzvcf AFR.vcf.gz --missing-site --out AFR
vcftools --gzvcf EUR.vcf.gz --missing-site --out EUR
vcftools --gzvcf AMR.vcf.gz --missing-site --out AMR
vcftools --gzvcf SAS.vcf.gz --missing-site --out SAS
vcftools --gzvcf EAS.vcf.gz --missing-site --out EAS
vcftools --gzvcf OCE.vcf.gz --missing-site --out OCE
vcftools --gzvcf MDE.vcf.gz --missing-site --out MDE
echo "..."
echo gathering freq stats
vcftools --gzvcf $isec --freq --out pooled
vcftools --gzvcf HGDP.vcf.gz --freq --out HGDP
vcftools --gzvcf ORIEN.vcf.gz --freq --out ORIEN
vcftools --gzvcf 1000GP.vcf.gz --freq --out 1000GP
vcftools --gzvcf TCGA.vcf.gz --freq --out TCGA
vcftools --gzvcf AFR.vcf.gz --freq --out AFR
vcftools --gzvcf EUR.vcf.gz --freq --out EUR
vcftools --gzvcf AMR.vcf.gz --freq --out AMR
vcftools --gzvcf SAS.vcf.gz --freq --out SAS
vcftools --gzvcf EAS.vcf.gz --freq --out EAS
vcftools --gzvcf OCE.vcf.gz --freq --out OCE
vcftools --gzvcf MDE.vcf.gz --freq --out MDE
echo "..."
echo gathering depth stats
vcftools --gzvcf $isec --site-mean-depth --out pooled
vcftools --gzvcf HGDP.vcf.gz --site-mean-depth --out HGDP
vcftools --gzvcf ORIEN.vcf.gz --site-mean-depth --out ORIEN
vcftools --gzvcf 1000GP.vcf.gz --site-mean-depth --out 1000GP
vcftools --gzvcf TCGA.vcf.gz --site-mean-depth --out TCGA
vcftools --gzvcf AFR.vcf.gz --site-mean-depth --out AFR
vcftools --gzvcf EUR.vcf.gz --site-mean-depth --out EUR
vcftools --gzvcf AMR.vcf.gz --site-mean-depth --out AMR
vcftools --gzvcf SAS.vcf.gz --site-mean-depth --out SAS
vcftools --gzvcf EAS.vcf.gz --site-mean-depth --out EAS
vcftools --gzvcf OCE.vcf.gz --site-mean-depth --out OCE
vcftools --gzvcf MDE.vcf.gz --site-mean-depth --out MDE
echo "..."
#echo gathering LD stats
#vcftools --gzvcf $isec --hap-r2 --out pooled
#vcftools--gzvcf HGDP.vcf.gz --hap-r2 --out HGDP
#vcftools --gzvcf ORIEN.vcf.gz --hap-r2 --out ORIEN
#vcftools --gzvcf 1000GP.vcf.gz --hap-r2 --out 1000GP
#vcftools --gzvcf TCGA.vcf.gz --hap-r2 --out TCGA
#vcftools --gzvcf AFR.vcf.gz --hap-r2 --out AFR
#vcftools --gzvcf EUR.vcf.gz --hap-r2 --out EUR
#vcftools --gzvcf AMR.vcf.gz --hap-r2 --out AMR
#vcftools --gzvcf SAS.vcf.gz --hap-r2 --out SAS
#vcftools --gzvcf EAS.vcf.gz --hap-r2 --out EAS
#vcftools --gzvcf OCE.vcf.gz --hap-r2 --out OCE
#vcftools --gzvcf MDE.vcf.gz --hap-r2 --out MDE
