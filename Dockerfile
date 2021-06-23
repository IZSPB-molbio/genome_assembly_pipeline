FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="10c12282e8190cfc4a29c2b1df72e0e9d0d0dd86c5d80949ca8df5719bd5536c"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/abricate.yml
#   prefix: /conda-envs/e98a0ca95f09435bd114cdb7154f5e9c
#   name: abricate
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - abricate
RUN mkdir -p /conda-envs/e98a0ca95f09435bd114cdb7154f5e9c
COPY envs/abricate.yml /conda-envs/e98a0ca95f09435bd114cdb7154f5e9c/environment.yaml

# Conda environment:
#   source: envs/prokka.yml
#   prefix: /conda-envs/03b1be2c3198629aa8ca18524890bc40
#   name: prokka
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - prokka
#     - blast=2.2
RUN mkdir -p /conda-envs/03b1be2c3198629aa8ca18524890bc40
COPY envs/prokka.yml /conda-envs/03b1be2c3198629aa8ca18524890bc40/environment.yaml

# Conda environment:
#   source: envs/quast.yml
#   prefix: /conda-envs/d0719b21ad0946113a5973cc51f2a10a
#   name: quast
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - quast
#     - bwa
RUN mkdir -p /conda-envs/d0719b21ad0946113a5973cc51f2a10a
COPY envs/quast.yml /conda-envs/d0719b21ad0946113a5973cc51f2a10a/environment.yaml

# Conda environment:
#   source: envs/r-stuff.yml
#   prefix: /conda-envs/64c37ec4bd79f97b48a13a1cc7ee66cf
#   name: r-stuff
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - r-tidyverse
#     - r-dt
RUN mkdir -p /conda-envs/64c37ec4bd79f97b48a13a1cc7ee66cf
COPY envs/r-stuff.yml /conda-envs/64c37ec4bd79f97b48a13a1cc7ee66cf/environment.yaml

# Conda environment:
#   source: envs/wgs.yml
#   prefix: /conda-envs/a8aed477430c1b6297f7fcaabfbed1dc
#   name: wgs
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - spades=3.15
#     - trimmomatic
#     - flash
#     - idba
#     - fastqc
#     - multiqc
#     - samtools
#     - bowtie2
#     - ivar
#     - pandas
#     - pilon
#   prefix: /lustrehome/dsimone/miniconda3/envs/wgs
RUN mkdir -p /conda-envs/a8aed477430c1b6297f7fcaabfbed1dc
COPY envs/wgs.yml /conda-envs/a8aed477430c1b6297f7fcaabfbed1dc/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/e98a0ca95f09435bd114cdb7154f5e9c --file /conda-envs/e98a0ca95f09435bd114cdb7154f5e9c/environment.yaml && \
    mamba env create --prefix /conda-envs/03b1be2c3198629aa8ca18524890bc40 --file /conda-envs/03b1be2c3198629aa8ca18524890bc40/environment.yaml && \
    mamba env create --prefix /conda-envs/d0719b21ad0946113a5973cc51f2a10a --file /conda-envs/d0719b21ad0946113a5973cc51f2a10a/environment.yaml && \
    mamba env create --prefix /conda-envs/64c37ec4bd79f97b48a13a1cc7ee66cf --file /conda-envs/64c37ec4bd79f97b48a13a1cc7ee66cf/environment.yaml && \
    mamba env create --prefix /conda-envs/a8aed477430c1b6297f7fcaabfbed1dc --file /conda-envs/a8aed477430c1b6297f7fcaabfbed1dc/environment.yaml && \
    mamba clean --all -y
