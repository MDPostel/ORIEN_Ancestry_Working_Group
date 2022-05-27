#!/usr/bin/bash
##SBATCH --partition=largemem
#SBATCH --partition=epyc-64
#SBATCH --mem=248G
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=bitriallelic_version

module load usc bzip2 bcftools vcftools

cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset
echo filtering ORIEN for biallelic and triallelic SNPs
bcftools view -m2 -M3 -v snps ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered.vcf.gz -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs.vcf.gz --threads $t

bcftools sort ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs.vcf.gz -T /scratch2/mpostel/ -m 490G -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted.vcf.gz --threads 64
bcftools view -e 'F_MISSING>=0.1' ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted.vcf.gz -Oz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10.vcf.gz --threads 64

#GET AN IDEA OF # OF SNPS PRIOR TO INTERSECTION W/ PLINK
module load plink
plink --vcf ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10.vcf.gz --make-bed --out ORIEN_missing10_bitriallelic

	#FILTER FOR BI+TRI-ALLELIC SNPS, MAF >=1%
echo filtering gnomAD for biallelic and triallelic SNPs and MAF >=1%
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/gnomAD/gnomad_hgdp_tgp_genomes_release3.1.2
for i in {1..22}; do bcftools view -m2 -M3 -v snps -e 'MAF<0.01' gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed.vcf.gz -Oz -o gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01.vcf.gz; bcftools index -t gnomad.genomes.v3.1.2.hgdp_tgp.chr${i}_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01.vcf.gz --threads $t; done


#7)---------------------------
#INTERSECT ORIEN W/ GNOMAD
cd /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset
echo running intersection
bcftools isec -p ./ ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10.vcf.gz /project/davidwcr_264/Projects/Mackenzie/ORIEN/gnomAD/gnomad_hgdp_tgp_genomes_release3.1.2/gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01.vcf.gz -Oz --threads 500
mv 0002.vcf.gz ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz --threads 64
mv 0003.vcf.gz gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01_isecORIEN.vcf.gz
bcftools index -t -f gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01_isecORIEN.vcf.gz --threads 64

plink --vcf ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz --make-bed --out ORIENisec_missing10_bitriallelic
plink --vcf gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01_isecORIEN.vcf.gz --make-bed --out GNOMADisec_bitriallelic --double-id

bcftools view -h ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz > ORIEN_header_fix_inbreeding.txt; sed -i 's/##INFO=<ID=InbreedingCoeff,Number=1/##INFO=<ID=InbreedingCoeff,Number=A/g' ORIEN_header_fix_inbreeding.txt

#Inbreeding coefficient definitions don't match in header; fix this
bcftools reheader -h ORIEN_header_fix_inbreeding.txt ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD.vcf.gz -o ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD_reheader.vcf.gz
bcftools index -t -f ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD_reheader.vcf.gz --threads 64
bcftools merge ORIEN_CombineGVCFs_Mar2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_sorted_missing10_isecGNOMAD_reheader.vcf.gz gnomad.genomes.v3.1.2.hgdp_tgp.chrALL_PASS_subsetbyextendedbed_bitriallelicSNPs_maf01_isecORIEN.vcf.gz -Oz -o gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_missing10_sorted_isec.vcf.gz --threads 64
bcftools index -t -f gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_missing10_sorted_isec.vcf.gz --threads 64

plink --vcf gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_missing10_sorted_isec.vcf.gz --make-bed --out ORIEN-GNOMAD-isec-missing10_bitriallelic --double-id
