#!/bin/bash

set -eo pipefail

sed -i 's/cpus = 8/cpus = 4/g' nextflow.config 
sed -i 's/cpus = 12/cpus = 4/g' nextflow.config
sed -i 's/cpus = 16/cpus = 4/g' nextflow.config 

export TERM=linux

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --fastq_input_long .github/data/fastq_long \
	 --db plassembler-db \
	 --outdir .github/data/test_output \
	 -with-report .github/data/test_output/nextflow_report.html \
 	 -with-trace .github/data/test_output/nextflow_trace.tsv
