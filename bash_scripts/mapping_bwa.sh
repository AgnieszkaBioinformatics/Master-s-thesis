#!/usr/bin/env bash
DIR="."
INDEX="./TAIR10_bwa_index/"
MAPPING="${DIR}/k27/mapped"
READS_DIR="${DIR}/k27/filtered_reads"
READS=("3h1_1_INPUT" "3H1_1" "3H1_2" "Ihp1_1_INPUT" "Ihp1_1" "Ihp1_3" "Mut_1_INPUT" "Mut_1" "Mut_2" "WT_1_INPUT" "WT_1" "WT_3")
THREADS=5
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bwa

for file in $READS; do
    bwa mem -P bwa_index "${READS_DIR}/${file}_R1_001.fastq" "${READS_DIR}/${file}_R2_001.fastq" > "${MAPPING}/${file}.sam"
done

