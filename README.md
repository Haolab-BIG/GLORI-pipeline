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

3.  **Download Files**:

      * `run_GLORI.sh`
      * `newGLORI.sif` (The Singularity container)
      * `get_anno/*.py`
      * `pipelines/*.py`
      * `run_GLORI.py`

4.  **Reference Data**:
        Generate annotation files (required).
        A directory containing bowtie/STAR index (Below are the detailed steps for the human hg38 genome. For other reference genomes, please download the corresponding files and replace them as needed).
      ```bash
      mkdir anno_data
      cd anno_data
      # Download Genome FASTA and GTF

      
