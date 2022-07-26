#!/usr/bin/bash

#singularity pull --name jdk.sif docker://openjdk
#singularity pull --name rfmix2.sif docker://indraniel/rfmix:v2
#singularity pull --name rfmix1.sif docker://indraniel/rfmix:v1

#|---------------------------------------RFMixver2---------------------------------------|

threads=64
singularity exec --bind ${dir}:/input ${dir}/jdk.sif java -jar beagle.25Mar22.4f6.jar gt=/input/{ORIEN_vcf} out={phased_ORIEN_vcf} nthreads=${threads}
singularity exec --bind ${dir}:/input ${dir}/jdk.sif java -jar beagle.25Mar22.4f6.jar gt=/input/{gnomAD_vcf} out={phased_gnomAD_vcf} nthreads=${threads}

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
singularity exec --bind ${dir}:/input ${dir}/rfmix2.sif RunRFMix.py -f ${query_vcf} -r ${ref_vcf} -m ${sample_map} -g ${genetic_map} -o ${out} --chromosome=${i} --n-threads=${threads} -G ${generations} --reanalyze-reference -e ${iterations} 
    done
    
#|---------------------------------------RFMixver1---------------------------------------|    

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
