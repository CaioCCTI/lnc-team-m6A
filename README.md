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
* [**Step 7: Run Minimap2 map to reference Pipeline on gene TUG1 and Run Samtools index**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-7-run-minimap2-map-to-reference-pipeline-on-gene-TUG1-and-run-samtools-index)
* [**Step 8: Run nanopolish eventalign generating their respective directories**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-8-run-nanopolish-eventalign-generating-their-respective-directories)
* [**Step 9: Run m6Anet dataprep step**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-9-run-m6anet-dataprep-step)
* [**Step 10: Run m6Anet Inference Step**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-10-run-m6anet-inference-step)
* [**Step 11: Visualize Results**](https://github.com/CaioCCTI/lnc-team-m6A?tab=readme-ov-file#step-11-visualize-results)

## Introduction

This project aims to function as a guide for scientists that are interested on working with analysis of lncRNAs (Long No Coding RNAs) and, in order to help with that, various bioinformatics tools (m6Anet oriented) and tutorials will be updated in this repository.

## Prerequisites
* [WSL](https://learn.microsoft.com/pt-br/windows/wsl/install) or [Ubuntu](https://ubuntu.com/download) machine
* [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent)

### Step 0: Find Long No Coding RNAs within literature

The main objective of this tutorial is to help find diferentially methylated DRACH motifs on long-no-coding RNAs related do SARS-CoV-2 infection. In order to do that, firstly, we recommend to search in the literature about genes that are LNC and seem to have significance on SARS-CoV-2 infection.

For further curation it is recommended to download and install [Geneious](https://www.geneious.com/) in order to use the "map-to-reference" tool of your reads against the gene of interest.

### Step 1: Obtain genes

In order to exemplify, firstly we will use the gene TUG1 available on: [https://www.ncbi.nlm.nih.gov/gene/823](https://www.ncbi.nlm.nih.gov/gene/823)

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

### Step 7: Run Minimap2 map to reference Pipeline on gene TUG1 and Run Samtools index
In this particular example, we are using Calu_48_Control to be mapped against TUG1 gene.

```sh
minimap2 -a -x map-ont TUG1.fasta ../../Calu_48_Control.fastq.gz | samtools view -bS -F 4 | samtools sort > TUG1.sorted.bam && samtools index TUG1.sorted.bam
```
### Step 8: Run nanopolish eventalign generating their respective directories

```sh
mkdir summary && mkdir eventalign && nanopolish eventalign --reads ../../Calu_48_Control.fastq.gz --bam TUG1.sorted.bam --genome TUG1.fasta --scale-events --signal-index --summary summary/summary.txt --threads 50 > eventalign/eventalign.txt
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

* "TUG1_control_data.site_proba.csv"
  
![image](https://github.com/user-attachments/assets/cf165a16-5095-4127-b01a-585dfa018e40)


* "TUG1_infected_data.site_proba.csv"

![image](https://github.com/user-attachments/assets/f637f19f-c028-4fd3-8ca2-31dd583f4050)

### Step 12: Analyzing Results

  In order to compare the distribuition of methylated sites per transcript in sequencing reads from uninfected and infected cells, we performed the Wilcoxon-Mann-Whitney (WMW) test, as implemented in R version 4.4.2 whithin the base R package. This nonparametric test is widely used to assess differences between two independent groups when their distributions do not meet the assumptions identical, whereas the alternative hypothesis (H₁) posits that the distributions differ significantly, specifically by detecting a difference in their medians.
  To apply the WMW test to our date, we firstly gathered all predictions generated by m6Anet and stored them in tables, which were subsequently imported into R as data frames and labeled as “uninfected” and “infected". Then, a probability theshold was applied (modification probability ≥ 0.6) to filter out sites with low confidence in m6A modification.
   The primary variable chosen to compare the methylation distribution between the two groups was the umber of reads per sites (n reads). Since the number of reads per site is not standardized by default, data normalization was required. For this, the number of reads per site in each group was divided by the total number of reads in that group using the equation expressed below. After normalization, the data underwent a negative logarithm transformation (−log) to enhance visualization consistency across samples:
   
![image](https://github.com/user-attachments/assets/2ec5bebf-453c-434e-8f35-d92ba4788f29) 

To visualize the distributional differences between groups, we generated violin plots using the ggplot2 package (version 3.5.1) in R. Violin plots provide comprehensive graphical representation of data distribution by combining a box plot with a kernel density estimate, enabling the visualization of both summary statistics and the underlying probability density function.
