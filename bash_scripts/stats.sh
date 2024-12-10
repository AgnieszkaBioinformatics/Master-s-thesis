#!/usr/bin/env
MAPPING="/home/am467109/mapping/v2/short"
source ~/miniconda3/etc/profile.d/conda.sh
THREADS=30

echo $(pwd)
for file in *sam
    do
        NAME=$(basename "$file" .sam)_stats
        #mkdir -p "$NAME"
        conda activate samtools
        echo $file
        samtools stats -d --threads $THREADS ${MAPPING}/${file} > ${MAPPING}/"${NAME}_result_s.txt"
        result=$(samtools view -@ $THREADS -c -q 60 ${file})
        echo "Result from samtools: $result"
        answer="${result} ${NAME}"
        echo "answer operation result: $answer"
        echo "$answer" >> "${MAPPING}/unique_long.txt"
    done

