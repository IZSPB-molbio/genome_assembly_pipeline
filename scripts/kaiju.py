import glob, os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

sample        = snakemake.wildcards.sample
R1            = snakemake.input.R1
R2            = snakemake.input.R2
outfile       = snakemake.output.kaiju_out
nodes_dmp     = snakemake.params.nodes_dmp
kaiju_db      = snakemake.params.kaiju_db
threads       = snakemake.threads
log           = snakemake.log[0]

with open(log, "w") as f:
    sys.stderr = sys.stdout = f
    
    shell("""
        kaiju \
        -t {nodes_dmp} \
        -f {kaiju_db} \
        -z {threads} \
        -o {outfile} \
        -i {R1} -j {R2}
        """)
