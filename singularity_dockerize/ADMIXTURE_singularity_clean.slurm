#!/usr/bin/bash

#-------------------------
  #ONLY INPUT USER CHANGES PER RUN:
#-------------------------
vcf=$(echo {ORIEN_vcf})
dir=$(echo $PWD)

#.........................................................................................................

#-------------------------
  #FILES THAT MUST BE IN WORKING DIRECTORY $dir:
#-------------------------
#1) SINGULARITY IMAGES:
   #singularity pull --name admixture.sif docker://evolbioinfo/admixture:v1.3.0
   #singularity pull --name vcf2bed.sif docker://dockerbiotools/bedops:latest
   #singularity pull bcftools.sif oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest
   #singularity pull --name plink.sif docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1
#2) gnomAD ADMIXTURE RESULTS (USE THIS PROJECTION EVERY TIME; DOES NOT CHANGE)
   #gnomAD=$(echo ORIEN.${K}.P.in)

#.........................................................................................................

#-------------------------
  #FYI: THIS IS HOW gnomAD INPUT WAS GENERATED (USE THIS PROJECTION EVERY TIME; DOES NOT CHANGE)
#-------------------------
#bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' gnomADonly_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_biallelicSNPs_missing10_sorted_isec.vcf.gz | bgzip -c > gnomAD.tsv.gz && tabix -s1 -b2 -e2 gnomAD.tsv.gz
#echo "Training Admixture K=${K} on reference individuals"
#refs=$(echo gnomAD.bed)
#singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture /input/${refs} ${K} -j${threads} --cv=10
#mv $HOME/gnomAD.${K}.P $dir/gnomAD.${K}.P; mv $HOME/gnomAD.${K}.Q $dir/gnomAD.${K}.Q
#cp gnomAD.${K}.P ORIEN.${K}.P.in

#.........................................................................................................

#-------------------------    
  #CONSTANTS:
#-------------------------
#K=20 #number of popuations; run separately for 5-20 or using for-loop as shown below
gnomAD=$(echo ORIEN.${K}.P.in)
threads=64
refs=$(echo gnomAD.tsv.gz)
#.........................................................................................................

#-------------------------
  #RUN ADMIXTURE:
#-------------------------
singularity exec --bind ${dir}:/input ${dir}/bcftools.sif bcftools view -T /input/${refs} /input/${vcf} -Oz -o /input/intersection.vcf.gz --threads 64
singularity exec --bind ${dir}:/input ${dir}/plink.sif plink1.9 --make-bed --vcf /input/intersection.vcf.gz --out /input/ORIEN
sample=$(echo ORIEN.bed)

export SINGULARITY_BINDPATH="$dir"

for K in {5..20}
    do
echo "Using K=${K} reference results for analysis of ORIEN samples"
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture -P /input/${sample} ${K} -j${threads} --cv=10
echo done
    done
