#!/usr/bin/env sh

# {sample} {genome_cds} {threads} {log}
sample="$1"
genome_cds="$2"
threads="$3"
log="$4"
output="$5" # results_dir/annotation/abricate/${sample}.done

base_abricate_outname=${output/.done/}

mkdir -p $(dirname $output)

for i in $(abricate --list | awk 'NR>1{print $1}')
do
    echo $i
    abricate \
    --db $i \
    --threads ${threads} \
    --minid 30 \
    --mincov 70 \
    ${genome_cds}  | awk -v genome_id=${sample} 'BEGIN{FS="\t";OFS="\t"}$1=genome_id{print $0}' > ${base_abricate_outname}_${i}.out
done 2> ${log} && touch ${output}
