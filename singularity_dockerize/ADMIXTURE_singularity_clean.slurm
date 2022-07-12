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
   #singularity pull bcftools.sif oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest
   #singularity pull --name plink.sif docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1
#2) gnomAD ADMIXTURE RESULTS (USE THIS PROJECTION EVERY TIME; DOES NOT CHANGE)
   #gnomAD=$(echo ORIEN.${K}.P.in)

#.........................................................................................................

#-------------------------
  #FYI: THIS IS HOW gnomAD INPUT WAS GENERATED (USE THIS PROJECTION EVERY TIME; DOES NOT CHANGE)
#-------------------------
#echo "Training Admixture K=${K} on reference individuals"
#refs=$(echo gnomAD.bed)
#singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture /input/${refs} ${K} -j${threads} --cv=10 -W ${dir}
#cp gnomAD.${K}.P ORIEN.${K}.P.in

#.........................................................................................................

#-------------------------    
  #CONSTANTS:
#-------------------------
#K=20 #number of popuations; run separately for 5-20 or using for-loop as shown below
gnomAD=$(echo ORIEN.${K}.P.in)
threads=64

#.........................................................................................................

#-------------------------
  #RUN ADMIXTURE:
#-------------------------
singularity exec --bind ${dir}:/input ${dir}/bcftools.sif bcftools isec /input/${refs} /input/${vcf} -Oz -o /input/intersection.vcf
singularity exec --bind ${dir}:/input ${dir}/plink.sif plink --make-bed --vcf /input/intersection.vcf --out /input/ORIEN.bed
sample=$(echo ORIEN.bed)

for K in {5..20}
    do
echo "Using K=${K} reference results for analysis of ORIEN samples"
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture -P /input/${sample} ${K} -j${threads} --cv=10 -W ${dir}
echo done
    done
