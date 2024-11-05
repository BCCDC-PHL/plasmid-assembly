#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys


def parse_flye_assembly_info(assembly_info_path):
    """
    Parse a Flye assembly info file into a list of dictionaries.

    :param assembly_info_path: path to the Flye assembly info file
    :type assembly_info_path: str
    :return: list of dictionaries containing assembly information
    :rtype: list[dict]
    """
    assembly_info = []
    field_name_translation = {
        'circ.': 'circular',
        'cov.': 'coverage',
        'mult.': 'multiplicity',
        'alt_group': 'alternative_group',
    }
    with open(assembly_info_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            assembly_info_row = {}
            for field, value in row.items():
                translated_field_name = field_name_translation.get(field, field)
                assembly_info_row[translated_field_name] = value
            try:
                assembly_info_row['length'] = int(assembly_info_row['length'])
            except ValueError as e:
                exit(f"Error: {e}")
            assembly_info.append(assembly_info_row)

    return assembly_info


def mark_chromosome_contigs(assembly_info, min_chromosome_contig_length):
    """
    Mark contigs as chromosomes based on length and circularity.

    :param assembly_info: List of dictionaries containing assembly information.
    :type assembly_info: list
    :param min_chromosome_contig_length: Minimum length of a contig to be considered a chromosome.
    :type min_chromosome_contig_length: int
    :return: List of dictionaries containing assembly information with a chromosome field.
    :rtype: list[dict]
    """
    for contig in assembly_info:
        if contig['length'] >= min_chromosome_contig_length and contig['circular'] == 'Y':
            contig['is_chromosome'] = True
        else:
            contig['is_chromosome'] = False

    return assembly_info


def parse_assembly_fasta(assembly_input_path):
    """
    Parse an assembly fasta file into a list of dictionaries.

    :param assembly_input_path: path to the assembly fasta file
    :type assembly_input_path: str
    :return: list of dictionaries containing assembly information
    :rtype: list[dict]
    """
    assembly = []
    with open(assembly_input_path, "r") as f:
        contig = {}
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if contig:
                    assembly.append(contig)
                contig = {}
                contig['seq_name'] = line.lstrip(">")
            else:
                if 'seq' in contig:
                    contig['seq'] += line
                else:
                    contig['seq'] = line

        if contig:
            assembly.append(contig)
                

        for contig in assembly:
            contig['length'] = len(contig['seq'])

    return assembly


def main(args):
    assembly_info = parse_flye_assembly_info(args.assembly_info)
    assembly_info = mark_chromosome_contigs(assembly_info, args.min_chromosome_contig_length)
    assembly = parse_assembly_fasta(args.assembly_input)

    for contig in assembly:
        seq_name = contig['seq_name']
        for contig_info in assembly_info:
            if contig['seq_name'] == contig_info['seq_name']:
                contig['circular'] = contig_info['circular']

    with open(args.output_chromosome_fasta, "w") as f:
        for contig in assembly:
            if contig['length'] >= args.min_chromosome_contig_length and contig['circular'] == 'Y':
                f.write(f">{contig['seq_name']} length={contig['length']} circular={contig['circular']}\n")
                for i in range(0, len(contig['seq']), 80):
                    f.write(f"{contig['seq'][i:i+80]}\n")

    if args.sample_id:
        for contig in assembly_info:
            contig['sample_id'] = args.sample_id
    output_fieldnames = [
        'sample_id',
        'seq_name',
        'length',
        'coverage',
        'circular',
        'is_chromosome',
        'repeat',
        'multiplicity',
        'alternative_group',
        'graph_path',
    ]
    with open(args.output_assembly_info, "w") as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames, extrasaction='ignore', dialect='unix', quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for contig in assembly_info:
            writer.writerow(contig)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", help="Sample ID")
    parser.add_argument("--assembly-info", help="Assembly info file from Flye")
    parser.add_argument("--assembly-input", help="Assembly fasta file")
    parser.add_argument("--min-chromosome-contig-length", type=int, default=1000000, help="Minimum length of a contig to be considered a chromosome (default: 1000000)")
    parser.add_argument("--output-chromosome-fasta", default='./chromosome.fa', help="Output fasta file containing chromosome contigs, (default: chromosome.fa)")
    parser.add_argument("--output-assembly-info", default='./assembly_info.csv', help="Output csv file containing assembly information, (default: assembly_info.csv)")
    args = parser.parse_args()
    main(args)
