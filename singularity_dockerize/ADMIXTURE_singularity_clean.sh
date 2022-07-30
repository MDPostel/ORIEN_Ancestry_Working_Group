#!/usr/bin/bash

#-------------------------
  #ONLY INPUT USER CHANGES PER RUN:
#-------------------------
samples=$(echo {ORIEN_vcf})
dir=$(echo $PWD)

#-------------------------
  #RUN ADMIXTURE:
#-------------------------
  #singularity pull --name admixture.sif docker://evolbioinfo/admixture:v1.3.0
  #singularity pull bcftools.sif oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest
  #singularity pull --name plink.sif docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1

threads=64
refs=$(echo gnomAD.vcf.gz)
export SINGULARITY_BINDPATH="$dir"
   
singularity exec --bind ${dir}:/input ${dir}/bcftools.sif bcftools isec /input/${refs}  /input/${vcf} -Oz  --threads 64 -p ./
mv 0002.vcf.gz refs.vcf.gz
mv 0003.vcf.gz samples.vcf.gz
singularity exec --bind ${dir}:/input ${dir}/bcftools.sif bcftools index -t -f /input/refs.vcf.gz --threads 64
singularity exec --bind ${dir}:/input ${dir}/bcftools.sif bcftools index -t -f /input/samples.vcf.gz --threads 64
singularity exec --bind ${dir}:/input ${dir}/plink.sif plink1.9 --make-bed --vcf /input/refs.vcf.gz --out /input/refs
singularity exec --bind ${dir}:/input ${dir}/plink.sif plink1.9 --make-bed --vcf /input/samples.vcf.gz --out /input/samples

#where refs = gnomAD, samples = ORIEN
for K in {5..20}
    do
echo "Training Admixture K=${K} on reference individuals"
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture /input/refs.vcf.gz ${K} -j${threads} --cv=10
cp refs.${K}.P samples.${K}.P.in
    done

for K in {5..20}
    do
echo "Using K=${K} reference results for analysis of ORIEN samples"
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture -P /input/samples.vcf.gz ${K} -j${threads} --cv=10
echo done
    done
