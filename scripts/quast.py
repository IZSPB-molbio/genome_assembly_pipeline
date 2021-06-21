import os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

assembly = snakemake.input.assembly
input_R1 = snakemake.input.R1
input_R2 = snakemake.input.R2
input_U  = snakemake.input.U
outdir   = snakemake.params.outdir

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
    {pe1} \
    {pe2} \
    {me} \
    -t 20 --glimmer \
    {assembly}
    """)
