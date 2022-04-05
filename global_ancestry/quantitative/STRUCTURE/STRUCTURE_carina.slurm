#!/bin/bash
##SBATCH --partition=largemem
#SBATCH --partition=oneweek
#SBATCH --mem=200GB
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=structure


cd ~/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/test_structure

mem=$(echo 200)

ulimit -n 4001

#conda init
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda creaate -n structure
#conda activate structure
#conda install pip
#conda install pgdspider
#conda update pgdspiderr
#conda install structure
#conda install bzip2 bcftools htslib

#VCFs subsetted by AIMs bed/vcf already for laser:
samples=$(echo ~/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/ORIEN_TCGA_laser_samples.vcf.gz)
refpops=$(echo ~/project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/1000GP_HGDP_laser_refs.vcf.gz)

#bcftools merge $samples $refpops -Oz -o ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz; tabix -f ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz
vcf=$(echo ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf.gz)
cp $vcf ORIEN_TCGA_1000GP_HGDP_STRUCTURE2.vcf.gz
gunzip ORIEN_TCGA_1000GP_HGDP_STRUCTURE2.vcf.gz; mv ORIEN_TCGA_1000GP_HGDP_STRUCTURE2.vcf ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf
vcf2=$(echo ORIEN_TCGA_1000GP_HGDP_STRUCTURE.vcf)

prefix=$(echo "ORIEN_TCGA_1000GP_HGDP_STRUCTURE")

#~/bin/plink-1.90/plink --vcf $vcf --recode structure --out $prefix --double-id
plink --vcf $vcf --write-snplist --double-id

#NUMLOCI=$(grep -o "chr[0-9][0-9]*" $vcf | grep -v "#" | wc -l)
#NUMINDS=$(bcftools view -h $vcf | grep -v "##" | sed 's/#//g' | sed 's/CHROM//g' |
 sed 's/POS//g' | sed 's/ID//g' | sed 's/REF//g' | sed 's/ALT//g' | sed 's/QUAL//g'
 | sed 's/FILTER//g' | sed 's/INFO//g' | sed 's/FORMAT//g' | sed '/^$/d' | wc -w)

NUMLOCI=$(grep -o "[0-9]* variants" plink.log | sort -u | cut -d" " -f1)
NUMINDS=$(grep -o "[0-9]* people" plink.log | sort -u | cut -d" " -f1)

for i in $(bcftools view -h $samples | grep -v "##" | sed 's/#//g' | sed 's/CHROM//
g' | sed 's/POS//g' | sed 's/ID//g' | sed 's/REF//g' | sed 's/ALT//g' | sed 's/QUAL
//g' | sed 's/FILTER//g' | sed 's/INFO//g' | sed 's/FORMAT//g' | tr "\t" "\n" | sed
 '/^$/d'); do echo "$i 0"; done > pops_ref_boolean.txt
for i in $(bcftools view -h $refpops | grep -v "##" | sed 's/#//g' | sed 's/CHROM//
g' | sed 's/POS//g' | sed 's/ID//g' | sed 's/REF//g' | sed 's/ALT//g' | sed 's/QUAL
//g' | sed 's/FILTER//g' | sed 's/INFO//g' | sed 's/FORMAT//g' | tr "\t" "\n" | sed
 '/^$/d'); do echo "$i 1"; done >> pops_ref_boolean.txt

sort -k1 pops_ref_boolean.txt > pops_ref_boolean_sorted.txt
#sort -k1 ${prefix}.recode.strct_in > ${prefix}.recode_sorted.strct_in
#join ${prefix}.recode_sorted.strct_in pops_ref_boolean_sorted.txt > ${prefix}_refs
1_samples0.recode.strct_in
#sed -e '1s/$/POPFLAG/' ${prefix}_refs1_samples0.recode.strct_in

echo "
# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=pop_defs.txt
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
# Only input following regions (refSeqName:start:end, multiple regions: whitespace
separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is bel
ow:
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
# Specify which data type should be included in the STRUCTURE file  (STRUCTURE can
only analyze one data type per file):
STRUCTURE_WRITER_DATA_TYPE_QUESTION=
# Save more specific fastSTRUCTURE format?
#STRUCTURE_WRITER_FAST_FORMAT_QUESTION=yes" > vcf2structure.spid
STRUCTURE_WRITER_FAST_FORMAT_QUESTION=no" > vcf2structure.spid

. "/home1/mpostel/anaconda3/condabin/conda"
conda activate structure
java -Xmx${mem}g -jar /home1/mpostel/anaconda3/envs/structure/share/pgdspider-2.1.1
.5-1/PGDSpider2.jar -inputfile $vcf2 -inputformat VCF -outputfile ${prefix}.recode.
strct_in -outputformat STRUCTURE -spid vcf2structure.spid

grep -v "SNP_" ORIEN_TCGA_1000GP_HGDP_STRUCTURE.recode.strct_in  > ORIEN_TCGA_1000G
P_HGDP_STRUCTURE.recode.strct_in_noheader

sed -i 's/\t/ /g' ORIEN_TCGA_1000GP_HGDP_STRUCTURE.recode.strct_in_noheader

IFS=$'\n'

grep "SNP_[0-9][0-9]*" ORIEN_TCGA_1000GP_HGDP_STRUCTURE.recode.strct_in > ORIEN_TCG
A_1000GP_HGDP_STRUCTURE.recode.strct_in_ADDEDPOPDEFS

for line in $(cat ORIEN_TCGA_1000GP_HGDP_STRUCTURE.recode.strct_in_noheader); do na
me=$(printf '%s\n' "$line" | cut -d" " -f1); match=$(grep $name pop_defs.txt | awk
'{print $2,$3}'); data=$(echo $line | awk '{for (f=2; f<=NF; ++f) { if (f!=2) {prin
tf("%s",OFS);} printf("%s",$f)}; printf "\n" }'); echo $name $match $data; done >>
ORIEN_TCGA_1000GP_HGDP_STRUCTURE.recode.strct_in_ADDEDPOPDEFS

#echo "
#define INFILE ${prefix}.recode.strct_in_ADDEDPOPDEFS
#define OUTFILE $prefix_out
#define NUMREPS 1000
#define MAXPOPS 5
#define BURNIN 10000
#define NUMLOCI 15783
#define USEPOPINFO 1
#define PFROMPOPFLAGONLY 1
#define NUMINDS 15772
#define LABEL 1
#define MARKERNAMES 1
#define POPDATA 1
#define POPFLAG 1
#define LOCDATA 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define ONEROWPERIND 0
#" > mainparams

#remove top line with periods from .strct_in?
#java -Xmx${mem}g -jar /home1/mpostel/bin/frontend/class/Structure.jar -K 5
#/home1/mpostel/bin/structure/console/structure -K 5

#FOR FASTSTRUCTURE
        #conda create -n structure python=2.7
        #conda activate structure
        #conda install structure
sed -i 's/[\t]*//' ${prefix}.recode.strct_in_ADDEDPOPDEFS
structure -m ./mainparams -K 5
