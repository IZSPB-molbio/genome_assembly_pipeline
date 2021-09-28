import glob, os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

sample        = snakemake.wildcards.sample
assembly_file = snakemake.input.assembly#.replace("pipeline_state/stage_9_terminate", "scaffolds.fasta")
outfile       = snakemake.output.res
#outdir        = snakemake.params.outdir
threads       = snakemake.threads
log           = snakemake.log[0]

with open(outfile, "w") as f:
    sys.stderr = sys.stdout = f
    
    shell("""
        referenceseeker \
        --threads {threads} \
        data/referenceseeker/bacteria-refseq \
        {assembly_file} 
        """)
