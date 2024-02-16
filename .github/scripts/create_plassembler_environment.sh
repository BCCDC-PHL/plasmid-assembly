#!/bin/bash

nextflow pull BCCDC-PHL/plasmid-assembly -r main

conda env create \
      -f ${HOME}/.nextflow/assets/BCCDC-PHL/plasmid-assembly/environments/plassembler.yml \
      -p ${HOME}/.conda/envs/plassembler-3ac96e6e6413c7c411c19f45d1796cea
