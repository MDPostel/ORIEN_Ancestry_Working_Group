#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=300G
#SBATCH --cpus-per-task=64
##SBATCH --ntasks=64
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=ORIEN_SNPs

ref=$(echo /project/davidwcr_392/ORIEN_Project/Resources/hs38DH.fa)
bed=$(echo "/project/davidwcr_392/ORIEN_Project/Resources/hg38_IDT_targets_sorted_merged_extended_sorted_merged_noalts_bedtoolsslop100bp.bed")
dbSNP=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/dbSNP_file/00-All_chr.vcf.gz")
gatk=$(echo /project/davidwcr_264/Packages/gatk/4.1.8.0/gatk)

mem=$(echo "290") #memory
t=$(echo "5") #threads
ulimit -n 4001 #files open
#old_date=$(echo August2021)
new_date=$(echo Mar2022)

module load usc openjdk python/3.9.2 bzip2 bcftools htslib

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset

#for each chromosome $i

#1)#---------------------------
/home1/mpostel/bin/gatk/gatk --java-options "-Xmx${mem}g -Xms${mem}g" CombineGVCFs \
   -L chr${i} \
   --tmp-dir /scratch1/mpostel \
   -R $ref \
   --variant {ORIEN_gvcf1.g.vcf.gz} \
# list all g.VCFs being used for SNP selection as shown above
   -O ORIEN_CombineGVCFs_Mar2022_chr${i}.vcf.gz


#2)#---------------------------
#GVCF <REF> BLOCKS to VCF FORMAT

	#OPTION 1: RUN GATK GENOTYPEGVCFS

echo "Genotyping g.VCFs"
name=$(echo ORIEN_CombineGVCFs_${new_date}_chr${i})

$gatk --java-options "-Xmx${mem}g -Xms${mem}g" GenotypeGVCFs \
   -R $ref \
   -V ORIEN_CombineGVCFs_${new_date}_chr${i}.vcf.gz \
   -O ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz \
   --verbosity ERROR \
   --allow-old-rms-mapping-quality-annotation-data
tabix -f ${name}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz

	#OPTION 2: RUN GATK COMBINEGVCFS
/home1/mpostel/bin/gatk/gatk --java-options "-Xmx${mem}g -Xms200g -XX:ParallelGCThreads=5" GenomicsDBImport \
   --genomicsdb-workspace-path /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/my_GenomicsDBImport_chr${i} \
   -L chr${i}  \
   --genomicsdb-vcf-buffer-size 50000 \
   --batch-size 50 \
   --reader-threads 5 \
   --tmp-dir /scratch2/mpostel \
   --sample-name-map orien_usc_plus_intermember_cohort.sample_map

#3)#---------------------------
#SUBSET BY EXOME BED AND DBSNP VCF; FILTER FOR DP>20

file=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/ORIEN_CombineGVCFs_${new_date}_chr${i}_GATKGenotypeGVCFs_gVCF2vcf.vcf.gz)
name=$(echo $file | cut -d"/" -f9 | cut -d"." -f1)

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/ORIEN

bcftools isec -R $bed -p ./ $file $dbSNP --threads $t -Oz
mv 0002.vcf.gz ${name}_subsetbed_dbSNPisec.vcf.gz
mkdir chr${i}; cd chr${i}
bcftools index -t -f ${name}_subsetbed_dbSNPisec.vcf.gz --threads $t

	#AS OF MAY 2022 - NO LONGER FILTERING ON DP
#ORIEN2=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/ORIEN/${name}_subsetbed_dbSNPisec.vcf.gz")
#bcftools filter --include 'INFO/DP>20' $ORIEN2 -Oz -o ${name}_subsetbed_dbSNPisec_DPgt20.vcf.gz --threads $t
#tabix -f ${name}_subsetbed_dbSNPisec_DPgt20.vcf.gz


#4)#---------------------------
#SORT

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/
file=$(echo "/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/isec/ORIEN/ORIEN_CombineGVCFs_${new_date}_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec.vcf.gz")

bcftools sort -m ${mem}G -T ./ $file -Oz -o ORIEN_CombineGVCFs_${new_date}_chr1_comb1_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_${new_date}_chr1_comb1_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_DPgt20_sorted.vcf.gz --threads $t


#5)---------------------------
#REMOVE UNDERSCORES FROM ORIEN SAMPLE NAMES -- WILL CAUSE DOWNSTREAM PROBLEMS IN PLINK OTHERWISE

bcftools view -h ORIEN_CombineGVCFs_${new_date}_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted.vcf.gz > ORIEN_header.txt
sed -i 's/_st_g//g' ORIEN_header.txt
bcftools reheader -h ORIEN_header.txt ORIEN_CombineGVCFs_${new_date}_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted.vcf.gz -o ORIEN_CombineGVCFs_${new_date}_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_${new_date}_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader.vcf.gz --threads $t


#6)---------------------------
#CONCAT CHROMOSOMES; FILTER FOR BIALLELIC SNPS AND MISSINGNESS <10%

for i in {1..22}
    do
echo "....."
echo "working on $i"
bcftools query -l ORIEN_CombineGVCFs_Mar2022_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted.vcf.gz | sort > samples_chr${i}.txt
bcftools view -S samples_chr${i}.txt ORIEN_CombineGVCFs_Mar2022_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted.vcf.gz > ORIEN_CombineGVCFs_Mar2022_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf
bgzip ORIEN_CombineGVCFs_Mar2022_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf.gz --threads 64
    done

for i in {1..22}; do echo ORIEN_CombineGVCFs_Mar2022_chr${i}_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf.gz; done > concat_list.txt

bcftools concat -f concat_list.txt -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf.gz --threads 64
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf.gz --threads 64
bcftools view -m2 -M2 -v snps ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf.gz -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs.vcf.gz --threads $t

bcftools sort ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs.vcf.gz -T /scratch2/mpostel/ -m 490G -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted.vcf.gz --threads 64
bcftools view -e 'F_MISSING>=0.1' ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted.vcf.gz -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10.vcf.gz --threads 64

#GET AN IDEA OF # OF SNPS PRIOR TO INTERSECTION W/ PLINK
module load plink
plink --vcf ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10.vcf.gz --make-bed --out ORIEN_missing10


#7)---------------------------
#INTERSECT ORIEN W/ GNOMAD

bcftools isec -p ./ ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10.vcf.gz /project/davidwcr_264/Projects/Mackenzie/ORIEN/gnomAD/gnomad_hgdp_tgp_genomes_release3.1.2/gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_snps_maf01.vcf.gz -Oz --threads 500
mv 0002.vcf.gz ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz --threads 64
mv 0003.vcf.gz gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_snps_maf01_isecORIEN.vcf.gz
bcftools index -t -f gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_snps_maf01_isecORIEN.vcf.gz --threads 64

plink --vcf ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz --make-bed --out ORIENisec_missing10
plink --vcf gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_snps_maf01_isecORIEN.vcf.gz --make-bed --out GNOMADisec --double-id

bcftools view -h ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz > ORIEN_header_fix_inbreeding.txt; sed -i 's/##INFO=<ID=InbreedingCoeff,Number=1/##INFO=<ID=InbreedingCoeff,Number=A/g' ORIEN_header_fix_inbreeding.txt

#Inbreeding coefficient definitions don't match in header; fix this
bcftools reheader -h ORIEN_header_fix_inbreeding.txt ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD_reheader.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD_reheader.vcf.gz --threads 64
bcftools merge ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_biallelicSNPs_sorted_missing10_isecGNOMAD_reheader.vcf.gz gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_snps_maf01_isecORIEN.vcf.gz -Oz -o gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_biallelicSNPs_missing10_sorted_isec.vcf.gz --threads 64
bcftools index -t -f gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_biallelicSNPs_missing10_sorted_isec.vcf.gz --threads 64

plink --vcf gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_biallelicSNPs_missing10_sorted_isec.vcf.gz --make-bed --out ORIEN-GNOMAD-isec-missing10 --double-id

