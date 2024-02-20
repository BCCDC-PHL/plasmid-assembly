#!/bin/bash


source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate plassembler-3ac96e6e6413c7c411c19f45d1796cea

plassembler download -d plassembler-db
