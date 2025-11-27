# GLORI-pipeline
This repository contains an integrated bioinformatics pipeline for analyzing GLORI (Global m6A Mapping with N3-Methyladenosine) sequencing data. Our pipeline is built upon and extends the tools provided by GLORI-tools(https://github.com/liucongcas/GLORI-tools).

GLORI provides absolute, single-base resolution quantification of m6A modifications. This automated workflow processes raw sequencing reads through quality control, alignment, and culminates in the precise identification of m6A sites.

# Part I Introduction
## i. Workflow
Here stands an throughout workflow of GLORI data analysis.
<img width="1194" height="288" alt="{493E7D28-5B20-47CB-A8C2-A25B7B26D575}" src="https://github.com/user-attachments/assets/000a6e54-2920-4c2a-a61a-412c2e599b75" />


## ii. Features
This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies. With a single command, the entire GLORI workflow—from raw FASTQ input through trimming, quality control, genome alignment, m6A sites calling —can be executed reproducibly on any compatible system.

# Part II Requirements
1.  **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```


4.  **Reference Data**:
        Generate annotation files (required).
        A directory containing bowtie/STAR index (Below are the detailed steps for the human hg38 genome. For other reference genomes, please download the corresponding files and replace them as needed).
      ```bash
      mkdir anno_files
      cd anno_files
      
      # 1. download files for annotation (required, using hg38 as example): 
      wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt 
      wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
      
      # 2. download reference genome and transcriptome
      wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
      wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz

      # Unzip the files
      gunzip GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
      gunzip hg38.fa.gz
      gunzip GCF_000001405.39_GRCh38.p13_rna.fna.gz
      ```

5.   **Required File Structure**

      ```bash
      project_directory/
      ├── Scripts
            ├── config.yaml
            └── GLORI.smk
      ├── Containers/
            └── newGLORI.sif
      ├── anno_files/
            ├── GCF_000001405.39_GRCh38.p13_genomic.gtf
            ├── GCF_000001405.39_GRCh38.p13_rna.fna
            └── hg38.fa
      ```
      
      - **GLORI.smk** — The main Snakemake workflow script.  
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.  
        ⚠️ Must be located in the same directory as `GLORI.smk`.
      - **newGLORI.sif** — Singularity container image with all required software and dependencies pre-installed.
      - **anno_files** — Reference genome and transcriptome; replace with your preferred reference.
      
# Part III Running

   * **Example code for ChIP-seq (Histomodification)**

      * **Step 1: Edit `config.yaml`**

        ```bash
       tool_dir: "/project_directory"  # Tool root directory
       input_dir: "/project_directory"  # Directory where raw Fastq files are located
       output_root: "/project_directory/output"  # Root directory for all output files
       sample_info: "/project_directory/sample_info.txt"  # Path to sample information file
       sif: "/project_directory/newGLORI.sif"  # Singularity image path
       threads: 20  # General number of threads (for downloading, index building, etc.)      
		   trim_galore_first:
           quality: 20        
           stringency: 1      
           error_rate: 0.3    
           min_length: 35  

       trim_galore_second:
           clip_r1: 10
           quality: 20
           min_length: 25
        ```

     * **Step 2: run snakemake**

        ```bash
        snakemake -s GLORI.smk \
        --use-singularity \
        --cores 8 \
        --singularity-args "--cleanenv --bind /project_directory:/project_directory" 
        ```

# Part IV Output

   * **Output Structure**
      ```bash
      output/
      ├──test_dir /
	     #SRR21356251_100k.fq.gz is test rawdata
	       ├──SRR21356251_100k_dup.fq.gz
         ├──SRR21356251_100k.fastq.gz_trimming_report.txt
         ├──SRR21356251_100k_rmdup.fq.gz
         ├──SRR21356251_100k_rmdup.fq.gz_trimming_report.txt
         ├──SRR21356251_100k_rmdup_trimmed.fastq
         ├──SRR21356251_100k_rmdup_trimmed_fastqc.html
         ├──SRR21356251_100k_rmdup_trimmed_fastqc.zip
         ├──SRR21356251_100k_rmdup_trimmed.fq.gz
         ├──SRR21356251_100k_trimmed_fastqc.html
         ├──SRR21356251_100k_trimmed_fastqc.zip
         └──SRR21356251_100k_trimmed.fq.gz
      ├──anno_files /
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno
         ├──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno.sorted
         └──GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base
      ├──new_annotated_Acites_output /
         ├──mapping-info/
         ├──testGLORI_prefix_A.bed_sorted
         ├──testGLORI_prefix_AGchanged_2.fq
         ├──testGLORI_prefix_merged.sorted.bam
         ├──testGLORI_prefix_merged.sorted.bam.bai
         ├──testGLORI_prefix.totalCR.txt
         └──testGLORI_prefix.totalformat.txt

        ```
    * **Output Interpretation**
      #new_annotated_Acites_output
      | Output files | Interpretation | 
      | :---: | :---: |
      | ${your_prefix}_merged.sorted.bam | The overall mapping of reads in .bam format |
      | ${your_prefix}.totalCR.txt | The text file containing the median value of the overall A-to-G conversion rate for each transcriptome and gene |
      | ${your_prefix}_A.bed_sorted | The m6A sites bed |


# Pipeline Reference

This pipeline is implemented with reference to the **GLORI-tools** methodology:
- **GLORI-tools**:https://github.com/liucongcas/GLORI-tools/
