#!/usr/bin/bash
#SBATCH --partition=largemem
#SBATCH --mem=300G
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=mpostel@usc.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-%j.out
#SBATCH --job-name=rfmix


threads=$SLURM_CPUS_PER_TASK
module load bzip2 bcftools htslib openjdk

for i in {1..22}
    do
chr=$(echo $i)
ulimit -n 4001

query_vcf=$(echo gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_missing10_sorted_isec_chr${chr}_ORIENonly_nochrnotation.vcf.gz)

ref_vcf=$(echo gnomad_ORIEN_May2022_chrALL_GATKGenotypeGVCFs_gVCF2vcf_subsetbed_dbSNPisec_ref-maf01_sorted_reheader_resorted_samplesordered_bitriallelicSNPs_missing10_sorted_isec_chr${chr}_gnomadonly_nochrnotation.vcf.gz)


. "/spack/apps/linux-centos7-x86_64/gcc-8.3.0/anaconda3-2019.10-bpb6unkhyvkirwkg44uqchcy5jyhzhvt/etc/profile.d/conda.sh"
##Phase the vcfs
#conda activate py-popgen
#vcf_phase.py --vcf $query_vcf --phase-algorithm beagle --out ORIEN_phased.vcf.gz --out-format vcf.gz
#vcf_phase.py --vcf $ref_vcf --phase-algorithm beagle --out gnomad_phased.vcf.gz --out-format vcf.gz
#tabix -f *phased.vcf.gz
#conda deactivate

conda activate rfmix

query_vcf2=$(echo ORIEN_phased.vcf.gz)
ref_vcf2=$(echo gnomad_phased.vcf.gz)
sample_map=$(echo refidpop2.map)

genetic_map=$(echo rfmix_hg38_chr${chr}.map)

out=$(echo rfmix_chr${chr})
generations=50 #defaults to 8
iterations=100 #default is off

#crf=9 ##to get script to finish faster; https://github.com/slowkoni/rfmix/issues/5
##use --reanalyze-reference if refs are also admixed

echo Running rfmix
rfmix -f ${query_vcf2} -r ${ref_vcf2} -m ${sample_map} -g ${genetic_map} -o ${out} --chromosome=${chr} --n-threads=${threads} -G ${generations} --reanalyze-reference -e ${iterations} #--crf-weight=${crf}

    done