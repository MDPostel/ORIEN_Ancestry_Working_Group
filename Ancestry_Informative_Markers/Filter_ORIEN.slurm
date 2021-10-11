#!/usr/bin/bash
#SBATCH --partition=dtg
#SBATCH --mem=60G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs_ORIEN

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ORIEN USC = SAMPLES
#---------------------------
#SET UP

ref=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa)
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged.bed")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
gatk=$(echo "/home1/mpostel/bin/gatk/gatk")

mem=$(echo "59") #memory
t=$(echo "600") #threads
ulimit -n 4001 #files open
date=$(echo August2021)

module load usc openjdk python/3.9.2 bzip2 bcftools htslib

#---------------------------
#GET THE SAMPLES (DOWNLOAD FROM DNANEXUS)
#samples saved in /project/davidwcr_392/ORIEN_Project/WES/g.VCFs
	# germline samples denoted with "st_g" (which is what we want). Don't use samples appended with "st_t" as these are TUMOR g.VCFs.
	
#---------------------------
#COMBINE GVCFS USING GATK

#for i in $(ls /project/davidwcr_392/ORIEN_Project/WES/g.VCFs/*st_g*gz); do echo "   --variant $i \\"; done
echo "Combining g.VCFs"
$gatk  --java-options "-Xmx${mem}g" CombineGVCFs  \
   -R $ref  \
   --tmp-dir /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset \
   --variant /project/davidwcr_392/ORIEN_Project/WES/g.VCFs/{sampleID}_st_g_markdup_recalibrated_Haplotyper.g.vcf.gz \
   # REPEAT ABOVE STEP FOR EACH SAMPLE ID
   -O ORIEN_CombineGVCFs_${date}.g.vcf.gz

tabix -f ORIEN_CombineGVCFs_${date}.g.vcf.gz

#---------------------------
#GVCF <REF> BLOCKS to VCF FORMAT

file=$(echo ORIEN_CombineGVCFs_${date}.g.vcf.gz)
echo "Genotyping merged g.VCF $file"
name=$(echo $file | cut -d"." -f1)
$gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \
   -R $ref \
   -V $file \
   -O ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz \
   --verbosity ERROR \
   --allow-old-rms-mapping-quality-annotation-data
#   -O ${name}_GATKGenotypeGVCFs.vcf.gz \
#tabix -f ${name}_GATKGenotypeGVCFs.vcf.gz
tabix -f ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz

#file=$(echo ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs.vcf.gz)
#echo "Expanding REF blocks into sites"
#name=$(echo $file | cut -d"." -f1)
#bcftools convert --threads $t --gvcf2vcf $file -f $ref -Oz -o ${name}_gVCF2vcf.vcf.gz
#tabix -f ${name}_gVCF2vcf.vcf.gz

#---------------------------
#SUBSET BY EXOME BED AND DBSNP VCF; FILTER FOR DP>20

file=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz)
name=$(echo $file | cut -d"/" -f9 | cut -d"." -f1)
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/ORIEN
bcftools isec -R $bed -p ./ $file $dbSNP --threads $t
mv 0002.vcf ${name}_subsetbed_dbSNPisec.vcf
bgzip ${name}_subsetbed_dbSNPisec.vcf
tabix -f ${name}_subsetbed_dbSNPisec.vcf.gz
ORIEN2=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/ORIEN/${name}_subsetbed_dbSNPisec.vcf.gz")
bcftools filter --include 'INFO/DP>20' $ORIEN2 -Oz -o ${name}_subsetbed_dbSNPisec_DPgt20.vcf.gz --threads $t
tabix -f ${name}_subsetbed_dbSNPisec_DPgt20.vcf.gz

#---------------------------
#SORT

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/
file=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/ORIEN/ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20.vcf.gz")

bcftools sort -m ${mem}G -T ./ $file -Oz -o ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz #--threads $t
tabix -f ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz

#---------------------------
#REMOVE UNDERSCORES FROM ORIEN SAMPLE NAMES -- WILL CAUSE DOWNSTREAM PROBLEMS IN PLINK OTHERWISE

bcftools view -h ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz > ORIEN_header.txt
sed -i 's/_st_g//g' ORIEN_header.txt
bcftools reheader -h ORIEN_header.txt ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz -o ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted_reheader.vcf.gz
tabix -f ORIEN_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted_reheader.vcf.gz

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
