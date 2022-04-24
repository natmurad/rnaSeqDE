FROM snakemake/snakemake:stable

# Install packages for RNA-seq Analysis
RUN conda install -c bioconda fastqc=0.11.9 trimmomatic=0.39 star=2.7.9a rsem=1.2.28

# Install multiqc manually
RUN git clone https://github.com/ewels/MultiQC.git && cd MultiQC && pip install .

# Download script folder from Github
RUN mkdir -p /home/rnaSeqPipe
RUN cd /home/rnaSeqPipe
RUN git clone https://github.com/natmurad/rnaSeqDE.git

# Download gdrive
RUN conda install -c conda-forge gdrive
