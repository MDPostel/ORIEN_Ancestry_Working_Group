#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=tcga


module load usc bzip2 bcftools htslib
t=$(echo 4000) #threads
mem=$(echo 199)

cd /project/davidwcr_666/TCGA/WES

ulimit -n 4001

#in case tool doesn't recognize .bai but does recognize .bam.bai:
  #for i in $(find /project/davidwcr_666/TCGA/WES/bam -type f -iname "*bai"); do rename -v bai bam.bai $i; done

#find /project/davidwcr_666/TCGA/WES/bam -type f -iname "*bam" > TCGA_germline_bamlist.txt
#find /project/davidwcr_264/Projects/Mackenzie/temp_TCGA -type f -iname "*bam" >> TCGA_germline_bamlist.txt
bamlist=$(echo "/project/davidwcr_666/TCGA/WES/TCGA_germline_bamlist.txt")
  #split up bamlists to parallelize

ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")
regions=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_merged_extended_merged.bed")

bcftools mpileup -b $bamlist1 -f $ref -R $regions --threads 4000 | bcftools call -m --threads 4000 -Oz -o TCGA_germline_AIMs_SNPs_bam2vcf_1.vcf.gz
tabix TCGA_germline_AIMs_SNPs_bam2vcf_1.vcf.gz

for i in {1..49}; do ls TCGA_germline_AIMs_SNPs_bam2vcf_smaller${i}.vcf.gz; done > TCGA_list.txt

#-----------------
#individual scripts

i=$(echo 1) #REPEAT 1-50 or however you want to break it up
bcftools view -R chr${i} TCGA_germline_AIMs_SNPs_bam2vcf_Sep2021.vcf.gz -Oz -o TCGA_germline_AIMs_SNPs_bam2vcf_Sep2021_chr${i}.g.vcf.gz
tabix -f TCGA_germline_AIMs_SNPs_bam2vcf_Sep2021_chr${i}.g.vcf.gz

#ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")

#---------------------------
#EXPAND REF BLOCKS

echo "expanding REF blocks"
file=$(echo TCGA_germline_AIMs_SNPs_bam2vcf_smaller${i}.vcf.gz)
name=$(echo $file | cut -d"." -f1)
bcftools convert $file --gvcf2vcf -f $ref -Oz -o ${name}_expandREFblocks.g.vcf.gz
tabix -f ${name}_expandREFblocks.g.vcf.gz

#---------------------------
#SUBSET BY EXOME BED AND DBSNP VCF; FILTER FOR DP>20

echo "subsetting by exome bed and dbSNP VCF; filter for DP>20"
file2=$(echo /project/davidwcr_666/TCGA/WES/${name}_expandREFblocks.g.vcf.gz)
name2=$(echo $file | cut -d"/" -f9 | cut -d"." -f1)
mkdir isec
cd isec
mkdir break_${i}; cd break_${i}
bcftools isec -R $bed $file2 $dbSNP -p ./ --threads $t
mv 0002.vcf ${name2}_subsetbed_dbSNPisec.vcf
bgzip ${name2}_subsetbed_dbSNPisec.vcf
tabix -f ${name2}_subsetbed_dbSNPisec.vcf.gz

TCGA2=$(echo "/project/davidwcr_666/TCGA/WES/isec/break_${i}/${name2}_subsetbed_dbSNPisec.vcf.gz")
bcftools filter --include 'INFO/DP>20' $TCGA2 -Oz -o ${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz --threads $t
tabix -f ${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz
                                                                                                                                                                                                   

#---------------------------
#CONCAT
cd /project/davidwcr_666/TCGA/WES

find /project/davidwcr_666/TCGA/WES -type f -iname *dbSNPisec_DPgt20.vcf.gz > TCGA_concat_list.txt

bcftools merge -l TCGA_concat_list.txt -g $ref -Oz -o TCGA_final.vcf.gz --threads $t
tabix -f TCGA_final.vcf.gz


#---------------------------
#SORT

bcftools sort -m ${mem}G -T ./ TCGA_final.vcf.gz -Oz -o TCGA_final_sorted_Oct.vcf.gz
tabix -f TCGA_final_sorted_Oct.vcf.gz
                                                                                                                                                                                                          
