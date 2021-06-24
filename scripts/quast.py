import os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

sample   = snakemake.wildcards.sample
assembly = snakemake.input.assembly.replace("pipeline_state/stage_9_terminate", "scaffolds.fasta")
input_R1 = snakemake.input.R1
input_R2 = snakemake.input.R2
input_U  = snakemake.input.U
outdir   = snakemake.params.outdir
threads  = snakemake.threads
log      = snakemake.log[0]

with open(log, "w") as f:
    sys.stderr = sys.stdout = f
    pe1 = ""
    pe2 = ""
    me  = ""
    
    for e in input_R1:
        pe1 += "--pe1 {e} ".format(e=e)
    
    for e in input_R2:
        pe2 += "--pe2 {e} ".format(e=e)
    
    for e in input_U:
        me  += "--single {e} ".format(e=e)
    
    shell("""
        quast \
        -o {outdir} \
        -l {sample} \
        {pe1} \
        {pe2} \
        {me} \
        -t {threads} --glimmer \
        {assembly}
        """)
