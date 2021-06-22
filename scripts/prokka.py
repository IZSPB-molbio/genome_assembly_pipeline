import glob, os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

sample        = snakemake.wildcards.sample
assembly_file = snakemake.input.assembly.replace("pipeline_state/stage_9_terminate", "scaffolds.fasta")
outdir        = snakemake.params.outdir
threads       = snakemake.threads

assembly_handle     = open(assembly_file, 'r')
new_assembly        = assembly_file.replace(".fasta", ".prokka.fasta")
new_assembly_handle = open(new_assembly, 'w')

for l in assembly_handle:
    if l.startswith(">"):
        l = "_".join(l.split("_")[:2]) + "\n"
    _ = new_assembly_handle.write(l)

new_assembly_handle.close()

shell("""
    prokka \
    --cpu {threads} \
    --outdir {outdir} \
    --prefix {sample} \
    {new_assembly}
    """)
