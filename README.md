# Welcome to LNC-Team repository!
Viral RNA-Related Analysis with a comprehensive pipeline introduction of our research workflow.

# Table of Contents  
* [**Introduction**](https://github.com/CaioCCTI/lnc-team/blob/main/README.md#introduction)

## Introduction

This project aims to function as a guide for scientists that are interested on working with analysis of lncRNAs (Long No Coding RNAs) and, in order to help with that, various bioinformatics tools and tutorials will be updated in this repository.

## Step 1: Obtain genes

In order to exemplify, firstly we will use the gene CAPN1 available on: [https://www.ncbi.nlm.nih.gov/gene/823](https://www.ncbi.nlm.nih.gov/gene/823)

## Step 2: Install m6Anet

In order to install m6Anet, you can follow the instructions on their repository ([GoekLab/m6Anet](https://github.com/GoekeLab/m6anet/blob/master/README.md?plain=1)) or execute one of the following commands under a python environment:

```sh
$ pip install m6anet
```
```sh
$ conda install m6anet
```
## Step 3: Install Nanopolish Eventalign

One of the prerequisites to run m6Anet is to install the Nanopolish Eventalign tool, since it will be important to index the fastq files that will be obtained shortly after.

```sh
git clone --recursive https://github.com/jts/nanopolish.git
cd nanopolish
make
```
