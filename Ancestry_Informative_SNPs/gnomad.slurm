#!/usr/bin/bash
#SBATCH --partition=epyc-64
##SBATCH --exclusive
#SBATCH --mem=240G
#SBATCH --cpus-per-task=64
##SBATCH --ntasks=64
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=gnomad

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/gnomAD/gnomad_hgdp_tgp_genomes_release3.1.2

#wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${#}.vcf.bgz

module load bzip2 bcftools htslib vcftools openjdk samtools jdk
t=64 #threads
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged_noalts_bedtoolsslop100bp.bed")

	#FILTER FOR PASS
echo filtering for PASS
for i in $(ls gnomad.genomes.v3.1.2.hgdp_tgp.chr*.vcf.gz); do name=$(echo $i | sed 's/.vcf.gz//g'); bcftools view -i '%FILTER=="PASS"' $i -Oz -o ${name}_PASS.vcf.gz --threads $t; bcftools index -t ${name}_PASS.vcf.gz --threads $t; done

	#SUBSET BY EXTENDED ORIEN BED
echo subsetting by extended ORIEN bed
for i in $(ls *PASS.vcf.gz); do prefix=$(echo $i | sed 's/.vcf.gz//g'); bcftools view -R $bed $i -Oz -o ${prefix}_subsetbyextendedbed.vcf.gz --threads $t; bcftools index -t -f ${prefix}_subsetbyextendedbed.vcf.gz --threads $t; done

	#FILTER FOR BIALLELIC SNPS, MAF >=1%
echo filtering for biallelic SNPs and MAF >=1%
for i in {1..22}; do bcftools view -m2 -M2 -v snps -e 'MAF<0.01' gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed.vcf.gz -Oz -o gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01.vcf.gz; bcftools index -t gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01.vcf.gz --threads $t; done

#AS OF MAY 2022 - NOT PERFORMING VQSR FILTER ON GNOMAD; ALREADY HAS QUALITY SCORE ANNOTATION "PASS"
	#VQSR
#echo VQSR
#GP1000_ref=$(echo "/project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa")
#GATK_ref=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.fasta)
#GP1000_ref2=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/GRCh38_full_analysis_set_plus_decoy_hla_fix.fa)
#HGDP_ref2=$(echo $GP1000_ref2)
#GP1000_phase3=$(echo ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz)
#dbSNP138=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/Homo_sapiens_assembly38.dbsnp138.vcf)
#dbSNP=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz)
#omni=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/1000G_omni2.5.hg38.vcf.gz)
#hapmap=$(echo /project/davidwcr_264/Projects/Mackenzie/resources/hapmap_3.3.hg38.vcf.gz)
#gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)

#for i in {1..22}
#    do
#echo running VariantAnnotator for chr${i}
#$gatk VariantAnnotator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
#   -V gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01.vcf.gz \
#   -O gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_QD.vcf.gz \
#   -A QualByDepth \
#   -A AS_QualByDepth \
#   -A RMSMappingQuality \
#   -A MappingQualityRankSumTest \
#   -A ReadPosRankSumTest \
#   -A FisherStrand \
#   -A StrandOddsRatio
#bcftools index -t -f gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_QD.vcf.gz --threads $t

#echo running VariantRecalibrator for chr${i}
#$gatk VariantRecalibrator --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
#   -R $GP1000_ref2 \
#   --variant gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_QD.vcf.gz \
#   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
#   --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
#   -an FS -an SOR \
#   -mode SNP \
#   --output gnomad_chr${i}.recal \
#   --tranches-file gnomad_chr${i}.tranches \
#   -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 98.5 -tranche 98.0 -tranche 97.5 -tranche 97.0 -tranche 96.5 -tranche 96.0 -tranche 95.5 -tranche 95.0 -tranche 94.5 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.5 -tranche 92.0 -tranche 91.5 -tranche 91.0 -tranche 90.5 -tranche 90.0 \
#   --rscript-file gnomad_chr${i}.plots.R

#echo running ApplyVQSR for chr${i}
#$gatk ApplyVQSR \
#   -R $GP1000_ref2 \
#   --variant gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_QD.vcf.gz \
#   --output gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_VQSR.vcf \
#   --tranches-file gnomad_chr${i}.tranches \
#   --truth-sensitivity-filter-level 90.0 \
#   --recal-file gnomad_chr${i}.recal \
#   -mode SNP

#sed -i 's/AD,Number=R/AD,Number=./g' GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf
#grep "AD,Number=" GP1000_subsetbed_isecdbSNP_DPgt20_sorted_snps_maf1_missing20_VQSR.vcf

#bgzip -f gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_VQSR.vcf
#bcftools index -t -f gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_VQSR.vcf.gz --threads $t

#echo filtering Tranches for chr${i}
#bcftools view -e '%FILTER=="VQSRTrancheSNP99.90to99.95" || %FILTER=="VQSRTrancheSNP99.95to100.00+" || %FILTER=="VQSRTrancheSNP99.95to100.00"'
#bcftools index -t -f gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_snps_maf01_VQSRPASSplusTopTranches.vcf.gz --threads $t
    done
