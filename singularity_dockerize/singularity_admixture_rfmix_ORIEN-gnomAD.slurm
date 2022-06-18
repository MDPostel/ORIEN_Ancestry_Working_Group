#!/usr/bin/bash


#singularity pull --name admixture.sif docker://evolbioinfo/admixture:v1.3.0
#singularity pull --name rfmix.sif docker://indraniel/rfmix:v2


#|---------------------------------------ADMIXTURE---------------------------------------|
refs=$(echo {gnomAD_bed})
samples=$(echo {ORIEN_bed})
pop=$(echo {admixture_supervised_key_pop})
dir=$(echo $PWD)

for K in {1..20}
    do
echo "Training Admixture K=${K} on reference individuals"
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture ${refs} ${K} -j${threads} --cv=10
echo "Using K=${K} reference results for analysis of ORIEN samples"
cp {gnomAD_bed}.${K}.P {ORIEN_bed}.${K}.P.in
singularity exec --bind ${dir}:/input ${dir}/admixture.sif admixture -P ${samples} ${K} -j${threads} --cv=10
echo done
    done
 

#|---------------------------------------RFMix---------------------------------------|
query_vcf=$(echo {phased_ORIEN_vcf})
ref_vcf=$(echo {phased_gnomAD_vcf})
sample_map=$(echo {ref_id_pop_map})
genetic_map=$(echo rfmix_hg38_chr${i}.map)
out=$(echo rfmix_chr${i})
generations=50 #defaults to 8
iterations=100 #default is off
threads=64
#optional: --crf-weight=${crf}

for i in {1..22}
    do
echo "Running RFMix on chr${i}"
singularity exec --bind ${dir}:/input ${dir}/rfmix.sif rfmix -f ${query_vcf} -r ${ref_vcf} -m ${sample_map} -g ${genetic_map} -o ${out} --chromosome=${i} --n-threads=${threads} -G ${generations} --reanalyze-reference -e ${iterations} 
    done