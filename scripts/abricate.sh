#!/usr/bin/env sh

# {sample} {genome_cds} {threads} {log}
sample="$1"
genome_cds="$2"
threads="$3"
log="$4"

mkdir -p annotation/abricate
for i in $(abricate --list | awk 'NR>1{print $1}')
do
    echo $i
    abricate \
    --db $i \
    --threads ${threads} \
    --minid 30 \
    --mincov 70 \
    ${genome_cds}  | awk -v genome_id=${sample} 'BEGIN{FS="\t";OFS="\t"}$1=genome_id{print $0}' > annotation/abricate/${sample}_${i}.out
done 2> ${log} && touch annotation/abricate/${sample}.done
