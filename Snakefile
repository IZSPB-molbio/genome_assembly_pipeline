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
import tempfile

# from Bio import SeqIO
# import numpy as np
import pandas as pd

from modules.config_parsers import (
    fastqc_outputs, get_bed_files, get_datasets_for_symlinks,
    get_fasta_files, get_genome_files, get_genome_single_vcf_files,
    get_genome_single_vcf_index_files, get_genome_vcf_files, get_kaiju_kronaplots, get_mt_genomes, get_mt_fasta,
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

checkm_data = config["checkm"]["data_dir"]
referenceseeker_data = config["referenceseeker"]["data_dir"]

results_dir = config["results_dir"]
qc                 = os.path.join(results_dir, config["qc"]["qc"])
qc_fastqc_raw      = os.path.join(config["qc"]["fastqc"]["raw"])
qc_fastqc_filtered = os.path.join(config["qc"]["fastqc"]["filtered"])
fastqc_folders = {"raw" : qc_fastqc_raw, "filtered" : qc_fastqc_filtered}
#
reads_raw_source    = config["reads_raw_source"]
raw_outpath         = os.path.join(results_dir, config["reads_raw"])
trimmomatic_outpath = os.path.join(config["read_processing"]["trimmomatic"]["outdir"])
#
assembly_spades_outpath = os.path.join(results_dir, config["assembly"]["spades"]["outdir"])

# parse run tabs
analysis_tab, datasets_tab = parse_config_tabs(analysis_tab_file=analysis_tab_file, datasets_tab_file=datasets_tab_file)
sample_list = list(analysis_tab["sample"])

rule all:
    input:
        # os.path.join(results_dir, "reports/multiqc_report.html"),
        # os.path.join(results_dir, "reports/abricate.html"),
        os.path.join(results_dir, "all/all_rmlst.out"),
        os.path.join(results_dir, "all/all_mlst.out"),
        # checkm output
        os.path.join(results_dir, "all", "all_checkm.out"),
        # copy of all scaffolds file, for easy retrieving
        # and checkm processing
        expand(os.path.join(results_dir, "all", "{sample}_scaffolds.fasta"), sample=sample_list),
        # kronaplot (read taxonomic classification)
        get_kaiju_kronaplots(datasets_tab=datasets_tab, results_dir=results_dir),
        # expand(os.path.join(results_dir, "{sample}/reports/{sample}_kaiju.out.html"), sample=sample_list),
        # expand(os.path.join(results_dir, "{sample}/qc/qualimap/qualimap_report.html"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/reports/multiqc_report.html"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/reports/{sample}_abricate.html"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/qc/referenceseeker/{sample}.tab"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/annotation/abricate/{sample}.done"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/annotation/prokka/{sample}/{sample}.sqn"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/assembly/spades/{sample}/scaffolds.fasta"), sample=sample_list),
        expand(os.path.join(results_dir, "{sample}/qc/quast/{sample}/report.pdf"), sample=sample_list),
        get_symlinks(datasets_tab, analysis_tab=analysis_tab, infolder=reads_raw_source,
                     outfolder=raw_outpath),
        # fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="raw", fastqc_folders = fastqc_folders),
        # fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="filtered", fastqc_folders = fastqc_folders),
        expand(fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="raw", fastqc_folders = fastqc_folders, results_dir = results_dir)),
        expand(fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="filtered", fastqc_folders = fastqc_folders, results_dir = results_dir)),

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
        cd {raw_outpath}
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
        cd {raw_outpath}
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
        html_report_R1 =  os.path.join(results_dir, "{sample}", qc_fastqc_raw, "{sample}_{library}.R1_fastqc.html"), #.R1_fastqc.html
        html_report_R2 =  os.path.join(results_dir, "{sample}", qc_fastqc_raw, "{sample}_{library}.R2_fastqc.html"), 
    params:
        outDir = os.path.join(results_dir, "{sample}", qc_fastqc_raw)
    threads:
        2
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of raw read files {input} with {version}, {wildcards}"
    log:
        os.path.join(results_dir, "{sample}/logs/qc/data/reads/raw/{sample}_{library}.log")
    conda: "envs/wgs.yml"
    shell:
        """
        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} {input} &> {log}
        """

rule kaiju:
    """Taxonomic classification of raw reads for QC purposes. Only allowed option is db folder"""
    input:
        R1 = os.path.join(raw_outpath, "{sample}_{library}.R1.fastq.gz"),
        R2 = os.path.join(raw_outpath, "{sample}_{library}.R2.fastq.gz"),
    output:
        kaiju_out = temp(os.path.join(results_dir, "{sample}/qc/kaiju/{sample}_{library}_kaiju.out"))
    params:
        kaiju_db  = config['read_processing']['kaiju']['db'],
        nodes_dmp = config['read_processing']['kaiju']['nodes_dmp']
    threads:
        config['read_processing']['kaiju']['threads']
    log:
        os.path.join(results_dir, "{sample}/logs/qc/{sample}_{library}_kaiju.log")
    conda:
        "envs/kaiju.yml"
    script:
        "scripts/kaiju.py"

rule kaiju2krona:
    input:
        rules.kaiju.output.kaiju_out
    output:
        kaiju_out_krona = os.path.join(results_dir, "{sample}/qc/kaiju/{sample}_{library}_kaiju.out.html")
    params:
        nodes_dmp   = config['read_processing']['kaiju']['nodes_dmp'],
        names_dmp   = config['read_processing']['kaiju']['names_dmp'],
        krona_input = lambda wildcards, output: output.kaiju_out_krona.replace(".html", ".krona")
    conda:
        "envs/kaiju.yml"
    shell:
        "kaiju2krona \
        -t {params.nodes_dmp} \
        -n {params.names_dmp} \
        -i {input} \
        -o {params.krona_input} && \
        ktImportText \
        -o {output.kaiju_out_krona} \
        {params.krona_input}"

rule trimmomatic:
    """ QCing and cleaning reads """
    params:
        java_cmd = config['read_processing']['trimmomatic']['java_cmd'],
        #jar_file = config['read_processing']['trimmomatic']['jar_file'],
        mem = config['read_processing']['trimmomatic']['java_vm_mem'],
        options = config['read_processing']['trimmomatic']['options'],
        processing_options = config['read_processing']['trimmomatic']['processing_options'],
        adapters = config['read_processing']['trimmomatic']['adapters'],
        out1P = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.R1.fastq.gz"),
        out2P = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.R2.fastq.gz"),
        out1U = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.1U.fastq.gz"),
        out2U = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.2U.fastq.gz"),
    input:
        R1 = os.path.join(raw_outpath, "{sample}_{library}.R1.fastq.gz"),
        R2 = os.path.join(raw_outpath, "{sample}_{library}.R2.fastq.gz"),
        # R1 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library, infolder=raw_outpath)[0],
        # R2 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library, infolder=raw_outpath)[1],
    output:
        out1P = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.R1.fastq.gz"),
        out2P = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.R2.fastq.gz"),
        out1U = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.U.fastq.gz"),
    threads:
        config['read_processing']['trimmomatic']['threads']
    # version:
    #     subprocess.check_output("trimmomatic -version", shell=True)
    # message:
    #     "Filtering read dataset {wildcards.sample}_{wildcards.library} with Trimmomatic. {wildcards}" # v{version}"
    log:
        os.path.join(results_dir, "{sample}/logs/data/reads/filtered/{sample}_{library}_trimmomatic.log")
    conda: "envs/wgs.yml"
    shell:
        """
        # export tap=$(which trimmomatic | sed "s/bin\/trimmomatic/share\/trimmomatic\/adapters\/{params.adapters}/g")
        export tap="data/all_adapters.fasta"
        
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
        html_report_R1 = os.path.join(results_dir, "{sample}", qc_fastqc_filtered, "{sample}_{library}.R1_fastqc.html"),
        html_report_R2 = os.path.join(results_dir, "{sample}", qc_fastqc_filtered, "{sample}_{library}.R2_fastqc.html"),
        html_report_U  = os.path.join(results_dir, "{sample}", qc_fastqc_filtered, "{sample}_{library}.U_fastqc.html"),
    params:
        outDir = os.path.join(results_dir, "{sample}", qc_fastqc_filtered)
    threads:
        3
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of filtered read files {input} with {version}"
    log:
        os.path.join(results_dir, "{sample}/logs/qc/data/reads/filtered/{sample}_{library}.log")
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
        R1 = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.notCombined_1.fastq.gz"),
        R2 = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.notCombined_2.fastq.gz"),
        U = os.path.join(results_dir, "{sample}", trimmomatic_outpath, "{sample}_{library}.extendedFrags.fastq.gz"),
    params:
        outdir = lambda wildcards, output: os.path.split(output.R1)[0]
    log:
        os.path.join(results_dir, "{sample}/logs/data/reads/filtered/flash/{sample}_{library}.log")
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
        # R1 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="R1", infolder=trimmomatic_outpath),
        # R2 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="R2", infolder=trimmomatic_outpath),
        # U  = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab, sample=wildcards.sample, mate="U", infolder=trimmomatic_outpath)
        R1 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab,
                                                    sample=wildcards.sample,
                                                    mate="R1",
                                                    infolder=os.path.join(results_dir, wildcards.sample, trimmomatic_outpath)),
        R2 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab,
                                                    sample=wildcards.sample,
                                                    mate="R2",
                                                    infolder=os.path.join(results_dir, wildcards.sample, trimmomatic_outpath)),
        U  = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab,
                                                    sample=wildcards.sample,
                                                    mate="U",
                                                    infolder=os.path.join(results_dir, wildcards.sample, trimmomatic_outpath))
    output:
        final_file = os.path.join(results_dir, "{sample}/assembly/spades/{sample}/scaffolds.fasta")
    params:
        outdir = lambda wildcards, output: output.final_file.replace("/scaffolds.fasta", "")
    conda:
        "envs/wgs.yml"
    threads:
        10
    log:
        os.path.join(results_dir, "{sample}/logs/assembly/spades/{sample}.log")
    script:
        "scripts/spades.py"

rule create_checkm_inputs:
    input:
        assembly = rules.assembly.output.final_file,
    output:
        assembly = os.path.join(results_dir, "all", "{sample}_scaffolds.fasta")
    shell:
        """
        cp {input.assembly} {output.assembly}
        """

rule checkm:
    input:
        assemblies = expand(os.path.join(results_dir, "all", "{sample}_scaffolds.fasta"), sample=sample_list)
    output:
        checkm_output = os.path.join(results_dir, "all", "all_checkm.out")
    params:
        indir           = os.path.join(results_dir, "all"),
        outdir          = os.path.join(results_dir, "all", "checkm"),
        checkm_data_dir = checkm_data,
        ext             = "fasta"
    conda:
        "envs/checkm.yml"
    threads:
        20
    log:
        os.path.join(results_dir, "logs/qc/checkm/checkm.log")
    shell:
        """
        export CHECKM_DATA_PATH={params.checkm_data_dir}
        checkm lineage_wf \
            --tab_table -t {threads} --pplacer_threads {threads} \
            -f {output.checkm_output} \
            -x {params.ext} \
            {params.indir} \
            {params.outdir}
        """

rule assembly_qc:
    input:
        assembly = rules.assembly.output.final_file,
        R1 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab,
                                                    sample=wildcards.sample,
                                                    mate="R1",
                                                    infolder=os.path.join(results_dir, wildcards.sample, trimmomatic_outpath)),
        R2 = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab,
                                                    sample=wildcards.sample,
                                                    mate="R2",
                                                    infolder=os.path.join(results_dir, wildcards.sample, trimmomatic_outpath)),
        U  = lambda wildcards: get_files_assembly(datasets_tab=datasets_tab,
                                                    sample=wildcards.sample,
                                                    mate="U",
                                                    infolder=os.path.join(results_dir, wildcards.sample, trimmomatic_outpath))
    output:
        pdf = os.path.join(results_dir, "{sample}/qc/quast/{sample}/report.pdf")
    params:
        outdir = lambda wildcards, output: output.pdf.replace("/report.pdf", "")
    conda:
        "envs/quast.yml"
    threads:
        10
    log:
        os.path.join(results_dir, "{sample}/logs/qc/quast/{sample}.log")
    script:
        "scripts/quast.py"

rule annotation:
    input:
        assembly = rules.assembly.output.final_file,
    output:
        sqn = os.path.join(results_dir, "{sample}/annotation/prokka/{sample}/{sample}.sqn")
    params:
        outdir = lambda wildcards, output: os.path.split(output.sqn)[0]
    conda:
        "envs/prokka.yml"
    threads:
        10
    log:
        os.path.join(results_dir, "{sample}/logs/annotation/prokka/{sample}.log")
    script:
        "scripts/prokka.py"

# rule post_assembly_alignment:
#     """Create input for qualimap"""
#     input:
#         R1 = 
#         R2 = 
#         assembly =
#     output:
#         sam = 
#     threads: 5
#     conda: "envs/wgs.yml"
#     log: 
#     script:
#         "scripts/post_assembly_alignment.py"
# 
# #os.path.join(results_dir, "{sample}/qc/qualimap/qualimap_report.html")
# rule qualimap:
#     input:
#         alignment      = rules.post_assembly_alignment2bam.output.bam,
#         annotation_sqn = rules.annotation.output.sqn,
#     output:
#         qualimap_report = os.path.join(results_dir, "{sample}/qc/qualimap/qualimap_report.html")
#     conda:
#         "envs/qualimap.yml"
#     log:
#         os.path.join(results_dir, "{sample}/logs/qc/qualimap/qualimap.log")
#     script:
#         "scripts/qualimap.py"

rule referenceseeker:
    input:
        assembly = rules.assembly.output.final_file,
    output:
        res = os.path.join(results_dir, "{sample}/qc/referenceseeker/{sample}.tab")
    conda:
        "envs/referenceseeker.yml"
    params:
        referenceseeker_data_dir = referenceseeker_data,
    threads:
        10
    log:
        os.path.join(results_dir, "{sample}/logs/qc/referenceseeker/{sample}.log")
    script:
        "scripts/referenceseeker.py"

rule abricate:
    input:
        sqn = rules.annotation.output.sqn
    output:
        done = os.path.join(results_dir, "{sample}/annotation/abricate/{sample}.done")
    conda:
        "envs/abricate.yml"
    threads:
        4
    log:
        os.path.join(results_dir, "{sample}/logs/annotation/abricate/{sample}.log")
    script:
        "scripts/abricate.py"

rule abricate_summary:
    """Collate abricate results for single sample"""
    input:
        os.path.join(results_dir, "{sample}/annotation/abricate/{sample}.done")
    output:
        os.path.join(results_dir, "{sample}/reports/{sample}_abricate.html")
    params:
        abricate_res_dir = lambda wildcards, input: os.path.split(input[0])[0]
    log:
        os.path.join(results_dir, "{sample}/logs/annotation/abricate/aggregate.log")
    conda:
        "envs/r-stuff.yml"
    script:
        "scripts/abricate_aggregate.R"

rule multiqc_all:
    """MultiQC report for single file"""
    input:
        quast_out = rules.assembly_qc.output.pdf,
        prokka_out = rules.annotation.output.sqn
        # expand(os.path.join(results_dir, "{sample}/annotation/prokka/{sample}/{sample}.sqn"), sample=sample_list),
    output:
        report = os.path.join(results_dir, "{sample}/reports/multiqc_report.html")
    conda:
        "envs/wgs.yml"
    params:
        input_dir = lambda wildcards: os.path.join(results_dir, wildcards.sample),
        multiqc_res_dir = lambda wildcards, output: os.path.split(output[0])[0],
        log_dir = os.path.join(results_dir, "logs")
        # assembly_dir    = os.path.join(results_dir, 
    log:
        os.path.join(results_dir, "{sample}/logs/qc/multiqc.log")
    shell:
        """
        multiqc -f \
        -o {params.multiqc_res_dir} \
        {params.input_dir} &> {log}
        """

rule mlst:
    input:
        assembly = rules.assembly.output.final_file,
    output:
        mlst = os.path.join(results_dir, "{sample}/sequence_typing/mlst/{sample}_mlst.out")
    params:
        outdir = lambda wildcards, output: os.path.split(output.mlst)[0]
    conda:
        "envs/mlst.yml"
    threads:
        1
    log:
        os.path.join(results_dir, "{sample}/logs/sequence_typing/mlst/{sample}.log")
    script:
        "scripts/mlst.py"

rule mlst_collate:
    input:
        expand(os.path.join(results_dir, "{sample}/sequence_typing/mlst/{sample}_mlst.out"), sample=sample_list)
    output:
        os.path.join(results_dir, "all/all_mlst.out")
    params:
        mlst_res_dir = lambda wildcards, input: os.path.split(input[0])[0]
    log:
        os.path.join(results_dir, "logs/sequence_typing/mlst/aggregate.log")
    # conda:
    #     "envs/r-stuff.yml"
    shell:
        """
        cat {input} > {output}
        """

rule rmlst_api:
    input:
        assembly = rules.assembly.output.final_file,
    output:
        rmlst_json     = os.path.join(results_dir, "{sample}/sequence_typing/rmlst/{sample}_rmlst.json"),
        rmlst_tab      = os.path.join(results_dir, "{sample}/sequence_typing/rmlst/{sample}_rmlst.tab")
    params:
        outdir = lambda wildcards, output: os.path.split(output.rmlst_json)[0]
    # conda:
    #     "envs/mlst.yml"
    threads:
        1
    log:
        "logs/sequence_typing/rmlst/{sample}.log"
    script:
        "scripts/rmlst.py"

rule rmlst_collate:
    input:
        expand(os.path.join(results_dir, "{sample}/sequence_typing/rmlst/{sample}_rmlst.tab"), sample=sample_list)
    output:
        os.path.join(results_dir, "all/all_rmlst.out")
    params:
        mlst_res_dir = lambda wildcards, input: os.path.split(input[0])[0]
    log:
        os.path.join(results_dir, "logs/sequence_typing/rmlst/aggregate.log")
    # conda:
    #     "envs/r-stuff.yml"
    shell:
        """
        cat {input} > {output}
        """
