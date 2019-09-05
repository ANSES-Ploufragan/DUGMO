# DUGMO
Detection of Unknown Genetically Modified Organisms

## Introduction

DUGMO is a bioinformatics pipeline for the detection of GMOs, including unknown GMOs, based on Illumina paired-end sequencing data. In first steps, coding sequences (CDS) are aligned through two successive BLASTN against the host pangenome with relevant tuned parameters to discriminate CDSs belonging to the wild genome (wgCDS) from potential GMO coding sequences (pgmCDS). Then, Bray-Curtis distances are calculated between the wgCDSs and the pgmCDSs based on the difference of genomic vocabulary between these two sets. Finally, two machine learning methods, namely random forest and generalized linear model, are carried out to target true GMO CDS(s). 

## Installation with conda 

Git glone
```
git clone https://github.com/jhurel/DUGMO.git
cd DUGMO
```
Create a conda environnement 
```
conda create -n DUGMO snakemake
conda activate DUGMO
```

## Params in configDUGMO.json

Params name | Description
------------|------------
`R1 ` | Input file Read 1 in fq or fq.gz format
`R2 ` | Input file Read 2 in fq or fq.gz format
`RefGenome` | Input file of reference genome in fasta format
`PangenomeCDS` | Input file of pangenome CDS in .fa.gz format
`PangenomeWhole` | Input file of pangenome whole genome in .fa.gz format
`outputFolder` | path of output folder for all results
`threads` | Number of threads to use (default: 10)

## Quick start

Launch the snakemake workflow
```
snakemake -rpk --use-conda
```

## Output files

Foldername | Description
-----------|------------
`Cleaning_pipeline` | Cleaning pipeline with Shovill, Trimmomatic, BWA, Bowtie2, Blast and Prokka
`DataBankGMO` | Databank files
`GMO` | GMO potential files
`HostGenome` | Host genome files
`MachineLearning` | Results at the final step
`Rmes` | Rmes output files
