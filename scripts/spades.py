import os, sys
from snakemake import shell

# input_R1 = sys.argv[1]
# input_R2 = sys.argv[2]
# input_U  = sys.argv[3]
# outdir   = sys.argv[4]

input_R1 = snakemake.input.R1
input_R2 = snakemake.input.R2
input_U  = snakemake.input.U
outdir   = snakemake.params.outdir

if len(input_R1) > 1:
    pe1 = ""
    for x, i in enumerate(input_R1):
        pe1 += "--pe{n}-1 {inputfile} ".format(n=x+1, inputfile=i)
    pe2 = ""
    for x, i in enumerate(input_R2):
        pe2 += "--pe{n}-2 {inputfile} ".format(n=x+1, inputfile=i)
    me  = ""
    for x, i in enumerate(input_U):
        me += "--pe{n}-m {inputfile} ".format(n=x+1, inputfile=i)
    shell("""
        spades.py \
        --isolate \
        -t 20 \
        -o {outdir} \
        {pe1} {pe2} {me}
        """)
else:
    shell("""
        spades.py \
        --isolate \
        -t 20 \
        -o {outdir} \
        -1 {input_R1} \
        -2 {input_R2} \
        --merged {input_U}
        """)
