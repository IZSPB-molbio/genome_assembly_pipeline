import glob, os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

sample        = snakemake.wildcards.sample
assembly_file = snakemake.input.assembly#.replace("pipeline_state/stage_9_terminate", "scaffolds.fasta")
outdir        = snakemake.params.outdir
threads       = snakemake.threads
log           = snakemake.log[0]

with open(log, "w") as f:
    sys.stderr = sys.stdout = f
    
    # assembly_handle     = open(assembly_file, 'r')
    # new_assembly        = assembly_file.replace(".fasta", ".prokka.fasta")
    # new_assembly_handle = open(new_assembly, 'w')
    # 
    # for l in assembly_handle:
    #     if l.startswith(">"):
    #         l = "_".join(l.split("_")[:2]) + "\n"
    #     _ = new_assembly_handle.write(l)
    # 
    # new_assembly_handle.close()
    
    shell("""
        mlst \
        --label {sample} \
        --nopath \
        --novel {outdir}/{sample}_new_alleles.fa \
        {assembly_file} >> {outdir}/{sample}_mlst.out
        """)
        # prokka \
        # --strain {sample} \
        # --force \
        # --cpu {threads} \
        # --outdir {outdir} \
        # --prefix {sample} \
        # {new_assembly}
        # """)
