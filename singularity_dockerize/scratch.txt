#!/usr/bin/bash


#singularity pull --name admixture.sif docker://evolbioinfo/admixture:v1.3.0
#singularity pull --name jdk.sif docker://openjdk
#singularity pull --name rfmix2.sif docker://indraniel/rfmix:v2
#singularity pull --name rfmix1.sif docker://indraniel/rfmix:v1
#singularity pull bcftools.sif oras://registry.forgemia.inra.fr/gafl/singularity/bcftools/bcftools:latest
#singularity pull --name plink.sif docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1

#|---------------------------------------ADMIXTURE---------------------------------------|

#only thing user changes:
vcf=$(echo {ORIEN_vcf})

#FYI: This is how the gnomAD Admixture input was generated
    #echo "Training Admixture K=${K} on reference individuals"
    #refs=$(echo gnomAD.bed)
    #singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture ${refs} ${K} -j${threads} --cv=10
    #cp gnomAD.${K}.P ORIEN.${K}.P.in

#constants:
pop=$(echo {admixture_supervised_key_pop})
dir=$(echo $PWD)
#K=20 #run separately for 5-20 or using for-loop as shown below
gnomAD=$(echo ORIEN.${K}.P.in)
threads=64

singularity exec --bind ${dir}:/input ${dir}/bcftools.sif bcftools isec $refs $vcf -Oz -o intersection.vcf
singularity exec --bind ${dir}:/input ${dir}/plink.sif plink --make-bed --vcf intersection.vcf --out ORIEN.bed
sample=$(echo ORIEN.bed)

for K in {5..20}
    do
echo "Using K=${K} reference results for analysis of ORIEN samples"
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture -P ${sample} ${K} -j${threads} --cv=10
echo done
    done
 
#|---------------------------------------RFMix2---------------------------------------|

threads=64
singularity exec --bind ${dir}:/input ${dir}/jdk.sif java -jar beagle.25Mar22.4f6.jar gt={ORIEN_vcf} out={phased_ORIEN_vcf} nthreads=${threads}
singularity exec --bind ${dir}:/input ${dir}/jdk.sif java -jar beagle.25Mar22.4f6.jar gt={gnomAD_vcf} out={phased_gnomAD_vcf} nthreads=${threads}

query_vcf=$(echo {phased_ORIEN_vcf})
ref_vcf=$(echo {phased_gnomAD_vcf})
sample_map=$(echo {ref_id_pop_map})
genetic_map=$(echo rfmix_hg38_chr${i}.map)
out=$(echo rfmix_chr${i})
generations=50 #defaults to 8
iterations=100 #default is off

#optional: --crf-weight=${crf}

for i in {1..22}
    do
echo "Running RFMix on chr${i}"
singularity exec --bind ${dir}:/input ${dir}/rfmix2.sif rfmix -f ${query_vcf} -r ${ref_vcf} -m ${sample_map} -g ${genetic_map} -o ${out} --chromosome=${i} --n-threads=${threads} -G ${generations} --reanalyze-reference -e ${iterations} 
    done
    
#|---------------------------------------RFMix1---------------------------------------|    

generations=50 #defaults to 8
iterations=100 #default is off
threads=64
window_size=$(echo "0.2")
out=$(echo rfmix_chr${i})
alleles=$(echo )
classes=$(echo )
locations=$(echo )

for i in {1..22}
    do
singularity exec --bind ${dir}:/input rfmix1.sif /opt/hall-lab/python-2.7.15/bin/RunRFMix.py -G ${generations} -e ${iterations} -w ${window_size} --num-threads ${threads} --use-reference-panels-in-EM --forward-backward -o ${out} ${alleles} ${classes} ${locations}
    done
