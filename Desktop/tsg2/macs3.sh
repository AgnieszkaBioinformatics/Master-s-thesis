macs3 callpeak -t H3K4me3_PA03.bam -c input_PA03.bam --broad -g hs --broad-cutoff 0.1 --outdir ~/macs3_dir/
macs3 callpeak -t H3K4me3_PA05.bam -c input_PA05.bam --broad -g hs --broad-cutoff 0.1 --outdir ~/macs3_dir/


# converting BED 6+3 to regular BED file
cut -f1-6 -d$'\t' 03.broadPeak > 03.bed
cut -f1-6 -d$'\t' 05.broadPeak > 05.bed