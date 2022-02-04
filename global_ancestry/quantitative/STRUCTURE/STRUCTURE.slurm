#!/usr/bin/bash
##SBATCH --partition=epyc-64
#SBATCH --partition=gpu
#SBATCH --mem=190G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=strucK8

K=$(echo 8)
module load usc bzip2 bcftools htslib
cd /project/davidwcr_788/ORIEN_Intermember_Project/ancestry_tests/structure
ulimit -n 4001

samples=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/ORIEN_TCGA.vcf.gz)
refpops=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/HGDP_1000GP.vcf.gz)

#bcftools merge $samples $refpops -Oz -o ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz; tabix -f ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz
vcf=$(echo /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/HGDP_1000GP_ORIEN_TCGA_snps_missing10_isec.vcf.gz)
#cp $vcf ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz
#gunzip ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz
vcf2=$(echo ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf)

prefix=$(echo "ORIEN_TCGA_1000GP_HGDP_STRUCTURE")

##~/bin/plink-1.90/plink --vcf $vcf --recode structure --out $prefix --double-id
#~/bin/plink-1.90/plink --vcf $vcf --write-snplist --double-id

#NUMLOCI=$(grep -o "chr[0-9][0-9]*" $vcf | grep -v "#" | wc -l)
#NUMINDS=$(bcftools view -h $vcf | grep -v "##" | sed 's/#//g' | sed 's/CHROM//g' | sed 's/POS//g' | sed 's/ID//g' | sed 's/REF//g' | sed 's/ALT//g' | sed 's/QUAL//g' | sed 's/FILTER//g' | sed 's/INFO//g' | sed 's/FORMAT//g' | sed '/^$/d' | wc -w)

NUMLOCI=$(grep -o "[0-9]* variants" /project/davidwcr_788/ORIEN_Intermember_Project/ancestry_tests/structure/plink.log | sort -u | cut -d" " -f1)
NUMINDS=$(grep -o "[0-9]* people" /project/davidwcr_788/ORIEN_Intermember_Project/ancestry_tests/structure/plink.log | sort -u | cut -d" " -f1)

#for i in $(bcftools view -h $samples | grep -v "##" | sed 's/#//g' | sed 's/CHROM//g' | sed 's/POS//g' | sed 's/ID//g' | sed 's/REF//g' | sed 's/ALT//g' | sed 's/QUAL//g' | sed 's/FILTER//g' | sed 's/INFO//g' | sed 's/FORMAT//g' | tr "\t" "\n" | sed '/^$/d'); do echo "$i 0"; done > pops_ref_boolean.txt
#for i in $(bcftools view -h $refpops | grep -v "##" | sed 's/#//g' | sed 's/CHROM//g' | sed 's/POS//g' | sed 's/ID//g' | sed 's/REF//g' | sed 's/ALT//g' | sed 's/QUAL//g' | sed 's/FILTER//g' | sed 's/INFO//g' | sed 's/FORMAT//g' | tr "\t" "\n" | sed '/^$/d'); do echo "$i 1"; done >> pops_ref_boolean.txt

#sort -k1 pops_ref_boolean.txt > pops_ref_boolean_sorted.txt
##sort -k1 ${prefix}.recode.strct_in > ${prefix}.recode_sorted.strct_in
##join ${prefix}.recode_sorted.strct_in pops_ref_boolean_sorted.txt > ${prefix}_refs1_samples0.recode.strct_in
##sed -e '1s/$/POPFLAG/' ${prefix}_refs1_samples0.recode.strct_in

echo "
# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=pops_ref_boolean.txt
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=diploid
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=yes
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=no
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=yes
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=no

# STRUCTURE Writer questions
WRITER_FORMAT=STRUCTURE

# Specify the locus/locus combination you want to write to the STRUCTURE file:
STRUCTURE_WRITER_LOCUS_COMBINATION_QUESTION=
# Do you want to include inter-marker distances?
STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=no
# Specify which data type should be included in the STRUCTURE file  (STRUCTURE can only analyze one data type per file):
STRUCTURE_WRITER_DATA_TYPE_QUESTION=
# Save more specific fastSTRUCTURE format?
STRUCTURE_WRITER_FAST_FORMAT_QUESTION=yes" > vcf2structure.spi

#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda create -n pgdspider
#conda activate pgdspider
#conda install pip
#conda install pgdspider
#conda update pgdspider

. "/spack/apps/linux-centos7-x86_64/gcc-8.3.0/anaconda3-2019.10-bpb6unkhyvkirwkg44uqchcy5jyhzhvt/etc/profile.d/conda.sh"
#conda activate pgdspider
#java -Xmx400g -jar /home1/mpostel/.conda/envs/pgdspider/share/pgdspider-2.1.1.5-1/PGDSpider2-cli.jar -inputfile $vcf2 -inputformat VCF -outputfile ${prefix}.recode.strct_in -outputformat STRUCTURE -spid vcf2structure.spi

cd K${K}
#num loci and indiv gotten from prev run
echo "
#define INFILE ../ORIEN_TCGA_1000GP_HGDP_STRUCTURE.recode.strct_in
#define OUTFILE $prefix
#define NUMREPS 1000
#define MAXPOPS $K
#define BURNIN 10000
#define NUMLOCI 23954
#define USEPOPINFO 0 #used to set to 1
#define PFROMPOPFLAGONLY 1
#define NUMINDS $NUMINDS
#define LABEL 1
#define MARKERNAMES 1
#define POPDATA 0
#define POPFLAG 1
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define ONEROWPERIND 0" > mainparams

#to estimate K
echo "
#define INFERALPHA 1
#define COMPUTEPROB 1" > extraparams

#remove top line with periods from .strct_in
conda activate structure
structure -K $K
