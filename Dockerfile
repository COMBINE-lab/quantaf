# DO NOT CHANGE
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

workdir /tmp/docker-build/work/

shell [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
env TZ='Etc/UTC'
env LANG='en_US.UTF-8'

arg DEBIAN_FRONTEND=noninteractive

run apt-get update && \
    apt-get install -y openjdk-17-jdk-headless

run apt-get update && \
    apt-get install -y libarchive13

# Install Mambaforge
run apt-get update --yes && \
    apt-get install --yes curl && \
    curl \
        --location \
        --fail \
        --remote-name \
        https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/conda -u && \
    rm Mambaforge-Linux-x86_64.sh

# Set conda PATH
env PATH=/opt/conda/bin:$PATH

run mamba install -y -c anaconda git
run mamba install -y -c bioconda salmon gffread bedtools alevin-fry pyroe
run mamba install -y -c conda-forge cxx-compiler time

run apt-get update && \
    apt-get install -y wget git ca-certificates bzip2 time

# Latch SDK
# DO NOT REMOVE

run /opt/conda/bin/pip install --upgrade latch==2.39.0.dev8
run mkdir /opt/latch

# Copy workflow data (use .dockerignore to skip files)
copy . /root/

# Nextflow funsies
copy .latch/bin/nextflow /root/nextflow
copy .latch/.nextflow /root/.nextflow
copy .latch/nf_entrypoint.py /root/nf_entrypoint.py

copy input_files/3M-february-2018.txt.gz /root/
copy input_files/737K-august-2016.txt /root/

# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root
