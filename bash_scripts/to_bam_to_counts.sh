#!/usr/bin/env
HOME="/home/am467109/"
MAPPING=${HOME}/"exp2"/"mapping"
source ~/miniconda3/etc/profile.d/conda.sh
THREADS=30
GTF=${HOME}/"Araport11_GTF_genes_transposons.current.gtf"
short=("brm1_31" "brm1_37" "brm5_43")

is_in() {
    LIST=$1
    VALUE=$2

    for x in $LIST; do
        if [ "$x" = "$VALUE" ]; then
            return 0
        else
            continue
        fi
    done
    return 1
}


for file in *sam
do
    echo "Processing file: $file"
    NAME=$(basename "$file" .sam)
    is_in "${short[*]}" "$NAME"
    if [[ $? != 0 ]]; then
        echo "Skipping $NAME"
        continue
    fi
    echo "Activating pipeline"
    conda activate samtools
    samtools view -bh "${NAME}.sam" > "${NAME}.bam"
    samtools sort "${NAME}.bam" > "${NAME}_sorted_s.bam"
    conda deactivate
    conda activate featurecounts
    featureCounts -p -O -T $THREADS -a ${GTF} -o "${HOME}/counts/v2/${NAME}_s.txt" "${NAME}_sorted_s.bam"
    conda deactivate
done


