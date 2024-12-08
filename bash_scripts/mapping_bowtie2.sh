#!/usr/bin/env bash
DIR="."
INDEX="~/mnt/chr8/data/am467109/TAIR10_index/"
MAPPING="${DIR}/k27/mapped"
READS_DIR="${DIR}/k27/filtered_reads"
READS=("3h1_1_INPUT" "3H1_1" "3H1_2" "Ihp1_1_INPUT" "Ihp1_1" "Ihp1_3" "Mut_1_INPUT" "Mut_1" "Mut_2" "WT_1_INPUT" "WT_1" "WT_3")
THREADS=5
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bowtie2

for file in $READS; do
    bowtie2 -x ${INDEX} -1 "${READS_DIR}/${file}_R1_001.fastq" -2 "${READS_DIR}/${file}_R2_001.fastq" > "${MAPPING}/${file}.sam"
done

