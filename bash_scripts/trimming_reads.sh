#!/usr/bin/env
THREADS=30
F="./exp2/raw_data"
FILTERED="./exp2/filtered_reads"
ADAPTERS="./miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters"

# reads
R1=""
R2=""

# trimmomatic
source ~/miniconda3/etc/profile.d/conda.sh
conda activate trimmomatic
trimmomatic PE -threads $THREADS -phred33 ${F}/${R1}.fastq ${F}/${R2}.fastq \
${FILTERED}/${R1}_trimmed.fastq ${FILTERED}/${R1}_untrimmed.fastq \
${FILTERED}/${R2}_trimmed.fastq ${FILTERED}/${R2}_untrimmed.fastq \
ILLUMINACLIP:${ADAPTERS}/TruSeq3-PE.fa:2:30:7 HEADCROP:10 TRAILING:10

#echo "Trimmomatic finished"

# fastqc
conda deactivate
conda activate fastqc

echo "starting fastqc"
fastqc ${FILTERED}/"${R1}_trimmed.fastq" ${FILTERED}/"${R2}_trimmed.fastq"  -o ${FILTERED} -t $THREADS
