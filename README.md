# Welcome to LNC-Team repository!
Viral RNA-Related Analysis with a comprehensive pipeline introduction of our research workflow with m6Anet.

# Table of Contents  
* [**Introduction**](https://github.com/CaioCCTI/lnc-team/blob/main/README.md#introduction)
* [**Prerequisites**](https://github.com/CaioCCTI/lnc-team/blob/main/README.md#prerequisites)
* [**Step1: Obtain Genes**](https://github.com/CaioCCTI/lnc-team-m6A/blob/main/README.md#step-1-obtain-genes)
* [**Step2: Install m6Anet**](https://github.com/CaioCCTI/lnc-team-m6A#step-2-install-m6anet)
* [**Step3: Install Nanopolish**](https://github.com/CaioCCTI/lnc-team-m6A#step-3-install-nanopolish)
* [**Step4: Install Samtools**](https://github.com/CaioCCTI/lnc-team-m6A#step-4-install-samtools)
* [**Step5: Install minimap2**](https://github.com/CaioCCTI/lnc-team-m6A#step-5-install-minimap2)
* [**Step 6: Download Brute Fast5 Data from NCBI**](https://github.com/CaioCCTI/lnc-team-m6A#step-6-download-brute-fast5-data-from-ncbi)
* [**Step 7: Run Minimap2 map to reference Pipeline on gene CAPN1 and Run Samtools index**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-7-run-minimap2-map-to-reference-pipeline-on-gene-capn1-and-run-samtools-index)
* [**Step 8: Run nanopolish eventalign generating their respective directories**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-8-run-nanopolish-eventalign-generating-their-respective-directories)
* [**Step 9: Run m6Anet dataprep step**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-9-run-m6anet-dataprep-step)
* [**Step 10: Run m6Anet Inference Step**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-10-run-m6anet-inference-step)
* [**Step 11: Visualize Results**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-11-visualize-results)

## Introduction

This project aims to function as a guide for scientists that are interested on working with analysis of lncRNAs (Long No Coding RNAs) and, in order to help with that, various bioinformatics tools (m6Anet oriented) and tutorials will be updated in this repository.

## Prerequisites
* [WSL](https://learn.microsoft.com/pt-br/windows/wsl/install) or [Ubuntu](https://ubuntu.com/download) machine
* [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent)

### Step 1: Obtain genes

In order to exemplify, firstly we will use the gene CAPN1 available on: [https://www.ncbi.nlm.nih.gov/gene/823](https://www.ncbi.nlm.nih.gov/gene/823)

### Step 2: Install m6Anet

In order to install m6Anet, you can follow the instructions on their repository ([GoekeLab/m6Anet](https://github.com/GoekeLab/m6anet/blob/master/README.md)) or execute one of the following commands under a python environment:

```sh
pip install m6anet
```
```sh
conda install m6anet
```
### Step 3: Install Nanopolish

One of the prerequisites to run m6Anet is to install the Nanopolish tool, since it will be important to index the fastq files that will be obtained shortly after.


Conda installation (Requires Conda managed environment)

```sh
conda install bioconda::nanopolish
```

Traditional installation (requires GCC Compiler)

```sh
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
make
```

### Step 4: Install Samtools

Samtools will be important to manage the alignment (.bam) files, solving their respective indexes.

```sh
conda install bioconda::samtools
```

### Step 5: Install minimap2

Minimap2 is needed to generate Binary Alignment Files (.bam).

```sh
conda install bioconda::minimap2
```

### Step 6: Download Brute Fast5 Data from NCBI

In order to start analyzing the lncRNAs, it is crucial to have the brute Fast5 files from one of the NCBI SRA projects, in this particular example we will be using the project PRJNA608224, with the following SRAs:

* [Nanopore direct RNA sequencing of SARS-CoV-2:Calu control 48hpi (SRR13089335)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR13089335&display=data-access) - Uninfected Calu-3 Cells

* [Nanopore direct RNA sequencing of SARS-CoV-2:Calu infected 48hpi (SRR13089334)](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR13089334&display=data-access) - Infected Calu-3 Cells

### Step 7: Run Minimap2 map to reference Pipeline on gene CAPN1 and Run Samtools index
In this particular example, we are using Calu_48_Control to be mapped against CAPN1 gene.

```sh
minimap2 -a -x map-ont CAPN1.fasta ../../Calu_48_Control.fastq.gz | samtools view -bS -F 4 | samtools sort > CAPN1.sorted.bam && samtools index CAPN1.sorted.bam
```
### Step 8: Run nanopolish eventalign generating their respective directories

```sh
mkdir summary && mkdir eventalign && nanopolish eventalign --reads ../../Calu_48_Control.fastq.gz --bam CAPN1.sorted.bam --genome CAPN1.fasta --scale-events --signal-index --summary summary/summary.txt --threads 50 > eventalign/eventalign.txt
```
### Step 9: Run m6Anet dataprep step

modify the "--n_processes" parameter according with your processor thread count

```sh
m6anet dataprep --eventalign eventalign/eventalign.txt --out_dir output/ --n_processes 12
```

### Step 10: Run m6Anet Inference Step

```sh
m6anet inference --input_dir output/ --out_dir output_csv/
```

### Step 11: Visualize Results

After Step 10, there will be two ".csv" files generated which will be useful for Further analysis.

Those files should be in the following format:

* "capn1_control_data.indiv_proba.csv"
  
![image](https://github.com/user-attachments/assets/4e824373-1cdb-4c83-b275-b650c803c9b9)

* "capn1_control_data.site_proba.csv"

![image](https://github.com/user-attachments/assets/5411a575-a294-46ae-920f-4d265e0a3228)

