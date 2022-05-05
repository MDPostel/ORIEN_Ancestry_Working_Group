#!/usr/bin/bash

module load usc bedtools2
bedtools slop -i hg38_IDT_targets_sorted_merged_extended_sorted_merged_noalts.bed -g GRCh38.v32.gtf -l 100 -r 100 > hg38_IDT_targets_sorted_merged_extended_sorted_merged_noalts_bedtoolsslop100bp.bed