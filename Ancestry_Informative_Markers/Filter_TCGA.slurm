#!/usr/bin/bash
#SBATCH --partition=oneweek
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=tcga


module load usc; module load bzip2; module load bcftools; module load htslib
t=$(echo 4000) #threads
mem=$(echo 199)

cd /project/davidwcr_666/TCGA/WES

ulimit -n 4001

#in case tool doesn't recognize .bai but does recognize .bam.bai:
#for i in $(find /project/davidwcr_666/TCGA/WES/bam -type f -iname "*bai"); do rename -v bai bam.bai $i; done

#find /project/davidwcr_666/TCGA/WES/bam -type f -iname "*bam" > TCGA_germline_bamlist.txt
#find /project/davidwcr_264/Projects/Mackenzie/temp_TCGA -type f -iname "*bam" >> TCGA_germline_bamlist.txt
bamlist=$(echo "/project/davidwcr_666/TCGA/WES/TCGA_germline_bamlist.txt")

ref=$(echo "/project/davidwcr_264/Resources/references/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.fa")

regions=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged.bed")

for i in {1..22}
    do
bcftools mpileup -b $bamlist -f $ref -r chr${i} -T $regions --threads $t -g 20 | bcftools call -m --threads $t -g 20 -Oz -o TCGA_germline_AIMs_SNPs_bam2vcf_chr${i}.vcf.gz
tabix -f TCGA_germline_AIMs_SNPs_bam2vcf_chr${i}.vcf.gz

#---------------------------
#EXPAND REF BLOCKS

echo "expanding REF blocks"
file=$(echo TCGA_germline_AIMs_SNPs_bam2vcf_chr${i}.vcf.gz)
name=$(echo $file | cut -d"." -f1)
bcftools convert $file --gvcf2vcf -f $ref -Oz -o ${name}_expandREFblocks.g.vcf.gz
tabix -f ${name}_expandREFblocks.g.vcf.gz

#---------------------------
#SUBSET BY EXOME BED AND DBSNP VCF; FILTER FOR DP>20

echo "subsetting by exome bed and dbSNP VCF; filter for DP>20"
file2=$(echo /project/davidwcr_666/TCGA/WES/${name}_expandREFblocks.g.vcf.gz)
name2=$(echo $file2 | cut -d"/" -f9 | cut -d"." -f1)
mkdir isec
cd isec
mkdir chr${i}; cd chr${i}
bcftools isec -R $bed $file2 $dbSNP -p ./ --threads $t
mv 0002.vcf ${name2}_subsetbed_dbSNPisec.vcf
bgzip ${name2}_subsetbed_dbSNPisec.vcf
tabix -f ${name2}_subsetbed_dbSNPisec.vcf.gz

TCGA2=$(echo "/project/davidwcr_666/TCGA/WES/isec/chr${i}/${name2}_subsetbed_dbSNPisec.vcf.gz")
bcftools filter --include 'INFO/DP>20' $TCGA2 -Oz -o ${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz --threads $t
tabix -f ${name2}_subsetbed_dbSNPisec_DPgt20.vcf.gz
    done
                                                                                                                                                                                                             
#---------------------------
#CONCAT

cd /project/davidwcr_666/TCGA/WES

find /project/davidwcr_666/TCGA/WES -type f -iname *dbSNPisec_DPgt20.vcf.gz > TCGA_concat_list.txt

bcftools merge -l TCGA_concat_list.txt -g $ref -Oz -o TCGA_final.vcf.gz --threads $t
tabix -f TCGA_final.vcf.gz

bcftools sort -m ${mem}G -T ./ TCGA_final.vcf.gz -Oz -o TCGA_final_sorted_Oct.vcf.gz
tabix -f TCGA_final_sorted_Oct.vcf.gz                                                                                                                                                                                                   
