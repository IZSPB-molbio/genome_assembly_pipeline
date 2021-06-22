import glob, os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

sample        = snakemake.wildcards.sample
genome_cds    = snakemake.input.sqn.replace(".sqn", ".ffn")
# assembly_file = snakemake.input.assembly.replace("pipeline_state/stage_9_terminate", "scaffolds.fasta")
# outdir        = snakemake.params.outdir
threads       = snakemake.threads
log           = snakemake.log[0]

import sys

# with open(log, "w") as f:
#     sys.stderr = sys.stdout = f
shell("""
    mkdir -p annotation/abricate
    for i in $(abricate --list | awk 'NR>1{print $1}')
    do
        echo $i
        abricate \
        --db $i \
        --threads {threads} \
        --minid 30 \
        --mincov 70 \
        {genome_cds}  | awk -v genome_id={sample} 'BEGIN{FS="\t";OFS="\t"}$1=genome_id{print $0}' > annotation/abricate/{sample}_${i}.out
    done 2> {log} && touch annotation/abricate/{sample}.done
    """)

# f.close()

# log <- file(snakemake@log[[1]], open="wt")
# sink(log)
# [rest of script]
