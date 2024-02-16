#!/bin/bash


source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate badread

mkdir -p .github/data/fastq_long

while IFS=',' read -r sample_id assembly; do
    badread simulate \
	--seed ${seed} \
	--reference ${assembly} \
	--length 50000,5000 \
	--quantity 15x \
	--junk_reads 1 \
	--random_reads 1 \
	--chimeras 1 \
	> .github/data/fastq_long/${sample_id}_RL.fastq

    gzip -f .github/data/fastq_long/${sample_id}_RL.fastq

done < .github/data/reads_to_simulate.csv

