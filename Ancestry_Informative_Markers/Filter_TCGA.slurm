#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=AIMs_TCGA1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TCGA = SAMPLES
#---------------------------
#SET UP

ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged.bed")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
gatk=$(echo "/home1/mpostel/bin/gatk/gatk")

mem=$(echo "199") #memory
t=$(echo "4000") #threads
ulimit -n 4001 #files open
date=$(echo Oct2021)

module load usc
module load openjdk
module load python/3.9.2
module load bzip2
module load bcftools
module load htslib

#---------------------------
#GET THE SAMPLES (DOWNLOAD FROM GDC)

#---------------------------
#HC GVCFS USING GATK
cd /project/davidwcr_666/TCGA/WES

echo "Making g.VCFs"
for i in $(cat TCGA_germline_bamlist.txt)
    do
name=$(echo $i | cut -d"/" -f8)

$gatk --java-options "-Xmx${mem}g" HaplotypeCaller \
   -R $ref \
   --tmp-dir /project/davidwcr_666/TCGA/WES \
   -I $i \
   -O ${name}_TCGA.g.vcf.gz \
   -ERC GVCF \
   -G StandardAnnotation \
   -G AS_StandardAnnotation

tabix -f ${name}_TCGA.g.vcf.gz
    done
    
echo "Combining g.VCFs"
$gatk --java-options "-Xmx${mem}g" CombineGVCFs  \
   -R $ref  \
   --tmp-dir /project/davidwcr_666/TCGA/WES \
   --variant {g.vcf.gz} \
#LIST ALL SAMPLES AS ABOVE, ONE PER LINE
   -O TCGA_CombineGVCFs_${date}.g.vcf.gz 
tabix -f TCGA_CombineGVCFs_${date}.g.vcf.gz

#---------------------------
#GVCF <REF> BLOCKS to VCF FORMAT

file=$(echo TCGA_CombineGVCFs_${date}.g.vcf.gz)
echo "Genotyping merged g.VCF $file"
name=$(echo $file | cut -d"." -f1)
$gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \
   -R $ref \
   -V $file \
   -O ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz \
   --verbosity ERROR \
   --allow-old-rms-mapping-quality-annotation-data

tabix -f ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz


#---------------------------
#SUBSET BY EXOME BED AND DBSNP VCF; FILTER FOR DP>20

file=$(echo /project/davidwcr_666/TCGA/WES/TCGA_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz)
name=$(echo $file | cut -d"/" -f9 | cut -d"." -f1)
cd isec
bcftools isec -R $bed -p ./ $file $dbSNP --threads $t
mv 0002.vcf ${name}_subsetbed_dbSNPisec.vcf
bgzip ${name}_subsetbed_dbSNPisec.vcf
tabix -f ${name}_subsetbed_dbSNPisec.vcf.gz

TCGA2=$(echo "/project/davidwcr_666/TCGA/WES/isec/${name}_subsetbed_dbSNPisec.vcf.gz")
bcftools filter --include 'INFO/DP>20' $TCGA2 -Oz -o ${name}_subsetbed_dbSNPisec_DPgt20.vcf.gz --threads $t
tabix -f ${name}_subsetbed_dbSNPisec_DPgt20.vcf.gz

#---------------------------
#SORT

cd /project/davidwcr_666/TCGA/WES
file=$(echo "/project/davidwcr_666/TCGA/WES/isec/TCGA_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20.vcf.gz")

bcftools sort -m ${mem}G -T ./ $file -Oz -o TCGA_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz #--threads $t
tabix -f TCGA_CombineGVCFs_${date}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
