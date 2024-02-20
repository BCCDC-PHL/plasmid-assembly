#!/bin/bash

mkdir -p .github/data/{ncbi_datasets,assemblies}

curl -o .github/data/ncbi_datasets/GCF_002968455.1.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_002968455.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,SEQUENCE_REPORT"

unzip .github/data/ncbi_datasets/GCF_002968455.1.zip -d .github/data/ncbi_datasets/GCF_002968455.1 && rm .github/data/ncbi_datasets/GCF_002968455.1.zip

cp .github/data/ncbi_datasets/GCF_002968455.1/ncbi_dataset/data/GCF_002968455.1/GCF_002968455.1_ASM296845v1_genomic.fna .github/data/assemblies/GCF_002968455.1.fa

rm -r .github/data/ncbi_datasets
