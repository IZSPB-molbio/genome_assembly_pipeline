import bz2
import gzip
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time

# from Bio import SeqIO
# import numpy as np
import pandas as pd

from modules.config_parsers import (
    fastqc_outputs, get_bed_files, get_datasets_for_symlinks,
    get_fasta_files, get_genome_files, get_genome_single_vcf_files,
    get_genome_single_vcf_index_files, get_genome_vcf_files, get_mt_genomes, get_mt_fasta,
    get_sample_bamfiles, get_symlinks, parse_config_tabs
)

from modules.general import (
    check_tmp_dir, gapped_fasta2contigs, get_files_assembly, get_seq_name, sam_to_fastq,
    sam_cov_handle2gapped_fasta, trimmomatic_input
)

### Parse values from config.yaml
configfile: "config.yaml"
#
analysis_tab_file = config["analysis_tab_file"]
datasets_tab_file = config["datasets_tab_file"]

results_dir = config["results_dir"]
qc_fastqc_raw      = os.path.join(results_dir, config["qc"]["fastqc"]["raw"])
qc_fastqc_filtered = os.path.join(results_dir, config["qc"]["fastqc"]["filtered"])
fastqc_folders = {"raw" : qc_fastqc_raw, "filtered" : qc_fastqc_filtered}
#
reads_raw_source    = config["reads_raw_source"]
raw_outpath         = os.path.join(results_dir, config["reads_raw"])
trimmomatic_outpath = os.path.join(results_dir, config["read_processing"]["trimmomatic"]["outdir"])
#
assembly_spades_outpath = os.path.join(results_dir, config["assembly"]["spades"]["outdir"])

# parse run tabs
analysis_tab, datasets_tab = parse_config_tabs(analysis_tab_file=analysis_tab_file, datasets_tab_file=datasets_tab_file)
sample_list = list(analysis_tab["sample"])

rule all:
    input:
        os.path.join(results_dir, "reports/multiqc_report.html"),
        os.path.join(results_dir, "reports/abricate.html"),
        expand(os.path.join(results_dir, "annotation/abricate/{sample}.done"), sample=sample_list),
        expand(os.path.join(results_dir, "annotation/prokka/{sample}/{sample}.sqn"), sample=sample_list),
        expand(os.path.join(results_dir, "assembly/spades/{sample}/scaffolds.fasta"), sample=sample_list),
        expand(os.path.join(results_dir, "qc/quast/{sample}/report.pdf"), sample=sample_list),
        get_symlinks(datasets_tab, analysis_tab=analysis_tab, infolder=reads_raw_source,
                     outfolder=raw_outpath),
        fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="raw", fastqc_folders = fastqc_folders),
        fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="filtered", fastqc_folders = fastqc_folders),

rule symlink_libraries:
    input:
        R1 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         infolder=reads_raw_source, outfolder=raw_outpath,
                                                         library=wildcards.library, d="R1"),
        R2 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         infolder=reads_raw_source, outfolder=raw_outpath,
                                                         library=wildcards.library, d="R2")
    output:
        R1 = raw_outpath + "/{sample}_{library}.R1.fastq.gz",
        R2 = raw_outpath + "/{sample}_{library}.R2.fastq.gz",
    shell:
        """
        cd {{raw_outpath}}
        ln -sf {input.R1} $(basename {output.R1})
        ln -sf {input.R2} $(basename {output.R2})
        """

rule symlink_libraries_uncompressed:
    input:
        R1 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         library=wildcards.library, d="R1"),
        R2 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         library=wildcards.library, d="R2")
    output:
        R1 = raw_outpath + "/{sample}_{library}.R1.fastq",
        R2 = raw_outpath + "/{sample}_{library}.R2.fastq",
    shell:
        """
        cd {{raw_outpath}}
        ln -sf {input.R1} $(basename {output.R1})
        ln -sf {input.R2} $(basename {output.R2})
        """

rule fastqc_raw:
    input:
        R1 = os.path.join(raw_outpath, "{sample}_{library}.R1.fastq.gz"),
        R2 = os.path.join(raw_outpath, "{sample}_{library}.R2.fastq.gz"),
        # R1 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library, infolder=raw_outpath)[0],
        # R2 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library, infolder=raw_outpath)[1],
    output:
        html_report_R1 =  os.path.join(qc_fastqc_raw, "{sample}_{library}.R1_fastqc.html"), #.R1_fastqc.html
        html_report_R2 =  os.path.join(qc_fastqc_raw, "{sample}_{library}.R2_fastqc.html"), 
    params:
        outDir = qc_fastqc_raw,
    threads:
        2
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of raw read files {input} with {version}, {wildcards}"
    log:
        os.path.join(results_dir, "logs/qc/data/reads/raw/{sample}_{library}.log")
    conda: "envs/wgs.yml"
    shell:
        """
        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} {input} &> {log}
        """

rule trimmomatic:
    """ QCing and cleaning reads """
    params:
        java_cmd = config['read_processing']['trimmomatic']['java_cmd'],
        #jar_file = config['read_processing']['trimmomatic']['jar_file'],
        mem = config['read_processing']['trimmomatic']['java_vm_mem'],
        options = config['read_processing']['trimmomatic']['options'],
        processing_options = config['read_processing']['trimmomatic']['processing_options'],
        out1P = os.path.join(trimmomatic_outpath, "{sample}_{library}.R1.fastq.gz"),
        out2P = os.path.join(trimmomatic_outpath, "{sample}_{library}.R2.fastq.gz"),
        out1U = os.path.join(trimmomatic_outpath, "{sample}_{library}.1U.fastq.gz"),
        out2U = os.path.join(trimmomatic_outpath, "{sample}_{library}.2U.fastq.gz"),
    input:
        R1 = os.path.join(raw_outpath, "{sample}_{library}.R1.fastq.gz"),
        R2 = os.path.join(raw_outpath, "{sample}_{library}.R2.fastq.gz"),
        # R1 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library, infolder=raw_outpath)[0],
        # R2 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library, infolder=raw_outpath)[1],
    output:
        out1P = os.path.join(trimmomatic_outpath, "{sample}_{library}.R1.fastq.gz"),
        out2P = os.path.join(trimmomatic_outpath, "{sample}_{library}.R2.fastq.gz"),
        out1U = os.path.join(trimmomatic_outpath, "{sample}_{library}.U.fastq.gz"),
    threads:
        config['read_processing']['trimmomatic']['threads']
    # version:
    #     subprocess.check_output("trimmomatic -version", shell=True)
    # message:
    #     "Filtering read dataset {wildcards.sample}_{wildcards.library} with Trimmomatic. {wildcards}" # v{version}"
    log:
        os.path.join(results_dir, "logs/data/reads/filtered/{sample}_{library}_trimmomatic.log")
    conda: "envs/wgs.yml"
    shell:
        """
        export tap=$(which trimmomatic | sed 's/bin\/trimmomatic/share\/trimmomatic\/adapters\/TruSeq3-PE.fa/g')
        
        trimmomatic PE {params.options} \
        -threads {threads} {input.R1} {input.R2} \
        {params.out1P} {params.out1U} {params.out2P} {params.out2U} \
        ILLUMINACLIP:$tap:2:30:10 {params.processing_options} &> {log}
        
        zcat {params.out1U} {params.out2U} | gzip > {output.out1U} && rm {params.out1U} {params.out2U}
        """
    # run:
    #     #trimmomatic_adapters_path = get_trimmomatic_adapters_path()
    #     shell("export tap=$(which trimmomatic | sed 's/bin\/trimmomatic/share\/trimmomatic\/adapters\/TruSeq3-PE.fa/g'); trimmomatic PE {params.options} -threads {threads} {input.R1} {input.R2} {params.out1P} {params.out1U} {params.out2P} {params.out2U} ILLUMINACLIP:$tap:2:30:10 {params.processing_options} &> {log}")
    #     shell("zcat {params.out1U} {params.out2U} | gzip > {output.out1U} && rm {params.out1U} {params.out2U}")

rule fastqc_filtered:
    input:
        out1P = rules.trimmomatic.output.out1P,
        out2P = rules.trimmomatic.output.out2P,
        out1U = rules.trimmomatic.output.out1U,
    output:
        html_report_R1 = os.path.join(qc_fastqc_filtered, "{sample}_{library}.R1_fastqc.html"),
        html_report_R2 = os.path.join(qc_fastqc_filtered, "{sample}_{library}.R2_fastqc.html"),
        html_report_U  = os.path.join(qc_fastqc_filtered, "{sample}_{library}.U_fastqc.html"),
    params:
        outDir = qc_fastqc_filtered
    threads:
        3
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of filtered read files {input} with {version}"
    log:
        os.path.join(results_dir, "logs/qc/data/reads/filtered/{sample}_{library}.log")
    conda: "envs/wgs.yml"
    shell:
        """

        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} {input} &> {log}

        """

rule merge_PE:
    input:
        R1 = rules.trimmomatic.output.out1P,
        R2 = rules.trimmomatic.output.out2P,
    output:
        R1 = os.path.join(trimmomatic_outpath, "{sample}_{library}.notCombined_1.fastq.gz"),
        R2 = os.path.join(trimmomatic_outpath, "{sample}_{library}.notCombined_2.fastq.gz"),
        U = os.path.join(trimmomatic_outpath, "{sample}_{library}.extendedFrags.fastq.gz"),
    params:
        outdir = lambda wildcards, output: os.path.split(output.R1)[0]
    log:
        os.path.join(results_dir, "logs/data/reads/filtered/flash/{sample}_{library}.log")
    conda: "envs/wgs.yml"
    shell:
        """
        flash \
        --allow-outies \
        -z \
        -t 5 \
        -d {params.outdir} \
        -o {wildcards.sample}_{wildcards.library} \
        {input.R1} \
        {input.R2} &> {log}
        """

rule assembly:
    # for the input, need to get two lists: R1 and R2.
    # if len(R1) > 1, use standard command line with -1, -2, --merged
    # otherwise need to use --pe<#>-1, --pe<#>-2, --pe<#>-m
    # at the moment, SE reads are left behind. 
    input:
        R1 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="R1", infolder=trimmomatic_outpath),
        R2 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="R2", infolder=trimmomatic_outpath),
        U  = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="U", infolder=trimmomatic_outpath)
    output:
        final_file = os.path.join(assembly_spades_outpath, "{sample}/scaffolds.fasta")
    params:
        outdir = lambda wildcards, output: output.final_file.replace("/scaffolds.fasta", "")
    conda:
        "envs/wgs.yml"
    threads:
        10
    log:
        os.path.join(results_dir, "logs/assembly/spades/{sample}.log")
    script:
        "scripts/spades.py"

rule assembly_qc:
    input:
        assembly = rules.assembly.output.final_file,
        R1       = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="R1", infolder=trimmomatic_outpath),
        R2       = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="R2", infolder=trimmomatic_outpath),
        U        = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="U", infolder=trimmomatic_outpath),
    output:
        pdf = os.path.join(results_dir, "qc/quast/{sample}/report.pdf")
    params:
        outdir = lambda wildcards, output: output.pdf.replace("/report.pdf", "")
    conda:
        "envs/quast.yml"
    threads:
        10
    log:
        os.path.join(results_dir, "logs/qc/quast/{sample}.log")
    script:
        "scripts/quast.py"

rule annotation:
    input:
        assembly = rules.assembly.output.final_file,
    output:
        sqn = os.path.join(results_dir, "annotation/prokka/{sample}/{sample}.sqn")
    params:
        outdir = lambda wildcards, output: os.path.split(output.sqn)[0]
    conda:
        "envs/prokka.yml"
    threads:
        10
    log:
        os.path.join(results_dir, "logs/annotation/prokka/{sample}.log")
    script:
        "scripts/prokka.py"

rule abricate:
    input:
        sqn = rules.annotation.output.sqn
    output:
        done = os.path.join(results_dir, "annotation/abricate/{sample}.done")
    conda:
        "envs/abricate.yml"
    threads:
        4
    log:
        os.path.join(results_dir, "logs/annotation/abricate/{sample}.log")
    script:
        "scripts/abricate.py"

rule abricate_summary:
    input:
        expand(os.path.join(results_dir, "annotation/abricate/{sample}.done"), sample=sample_list)
    output:
        os.path.join(results_dir, "reports/abricate.html")
    params:
        abricate_res_dir = lambda wildcards, input: os.path.split(input[0])[0]
    log:
        os.path.join(results_dir, "logs/annotation/abricate/aggregate.log")
    conda:
        "envs/r-stuff.yml"
    script:
        "scripts/abricate_aggregate.R"

rule multiqc_all:
    input:
        expand(os.path.join(results_dir, "qc/quast/{sample}/report.pdf"), sample=sample_list),
        expand(os.path.join(results_dir, "annotation/prokka/{sample}/{sample}.sqn"), sample=sample_list),
    output:
        report = os.path.join(results_dir, "reports/multiqc_report.html")
    conda:
        "envs/wgs.yml"
    params:
        multiqc_res_dir = lambda wildcards, output: os.path.split(output[0])[0]
        # assembly_dir    = os.path.join(results_dir, 
    log:
        os.path.join(results_dir, "logs/qc/multiqc.log")
    shell:
        """
        multiqc -f \
        -x {assembly_spades_outpath} \
        -x {raw_outpath} \
        -o {params.multiqc_res_dir} .
        """
