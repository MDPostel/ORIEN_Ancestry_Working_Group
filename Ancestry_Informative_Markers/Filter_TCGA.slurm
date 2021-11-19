#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mem=1000G
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

mem=$(echo "990") #memory
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
    
#---------------------------
#Combining GVCFs -- exons

module load openjdk
picardjar=$(echo "/project/davidwcr_264/Packages/picard-tools/picard-2.9.0/picard.jar")
input_bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged.bed")
#TCGA_ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")
#cat $TCGA_ref | grep -v "KI" | grep -v "GL" | grep -v "chrUn" | grep -v "random" | grep -v HLA | grep -v alt | grep -v chrEBV | grep -v chrM > /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa
ref=$(echo /project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome_fix.fa)
dic=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hs38DH.dict")
java -jar $picardjar BedToIntervalList \
      I=$input_bed \
      O=ORIEN.interval_list \
      SD=$dict
grep -v "chrUn" ORIEN.interval_list | grep -v "_random" | grep -v "alt" | grep -v "HLA" | grep -v chrX | grep -v chrY | grep -v chrM | grep -v chrEBV > ORIEN_autosomes.interval_list
intervals=$(echo "/scratch/mpostel/ORIEN_autosomes.interval_list")

echo "Combining g.VCFs"
interval_list=$(echo /scratch/mpostel/ORIEN_autosomes.interval_list)
/home1/mpostel/bin/gatk/gatk --java-options "-Xmx${mem}g -Xms${mem}g" GenomicsDBImport \
    --genomicsdb-workspace-path /project/davidwcr_666/TCGA/WES/my_GenomicsDBImport \
    --tmp-dir /project/davidwcr_666/TCGA/WES \
    --intervals $interval_list \
    -V #list each individual TCGA gVCF on separate lines prefaced by "-V" and followed by "\" (n=11,244)
    --reader-threads 600 \
    --batch-size 250
    
#---------------------------
#JOINT-CALLING; GENOTYPE MERGED GVCF

file=$(echo TCGA_CombineGVCFs_${date}.g.vcf.gz)
echo "Genotyping merged g.VCF $file"
name=$(echo $file | cut -d"." -f1)

$gatk --java-options "-Xmx${mem}g" GenotypeGVCFs \
   -R $ref \
   -V gendb:///project/davidwcr_666/TCGA/WES/my_GenomicsDBImport \
   -O ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz \
   --verbosity ERROR \
   -G StandardAnnotation -newQual \
   -G AS_StandardAnnotation -newQual 
   
#   --allow-old-rms-mapping-quality-annotation-data

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
