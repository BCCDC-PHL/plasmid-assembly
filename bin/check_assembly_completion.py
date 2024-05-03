#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys

def parse_flye_assembly_info(assembly_info_path):
    """
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


def check_completeness(assembly_info, min_chromosome_contig_length):
    """
    """
    complete = False
    for contig in assembly_info:
        if contig['length'] >= min_chromosome_contig_length:
            complete = True
            break

    return complete


def main(args):
    assembly_info = parse_flye_assembly_info(args.input)
    assembly_complete = check_completeness(assembly_info, args.min_chromosome_contig_length)

    longest_contig = max(assembly_info, key=lambda x: x['length'])
    longest_contig_id = longest_contig.get('seq_name', None)
    longest_contig_circular = True if longest_contig.get('circular', None) == 'Y' else False

    output = {
        'sample_id': args.sample_id,
        'assembly_completeness': 'complete' if assembly_complete else 'incomplete',
        'num_total_contigs': len(assembly_info),
        'num_circular_contigs': len([contig for contig in assembly_info if contig.get('circular', None) == 'Y']),
        'longest_contig_id': longest_contig_id,
        'longest_contig_length': longest_contig['length'],
        'longest_contig_circular': longest_contig_circular,
    }

    writer = csv.DictWriter(sys.stdout, fieldnames=output.keys(), dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    writer.writerow(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='assembly summary file from flye')
    parser.add_argument('-s', '--sample-id', help='sample id')
    parser.add_argument('--min-chromosome-contig-length', type=int, default=1000000, help='Minimum length of a contig to be considered a chromosome (default: 1000000)')
    args = parser.parse_args()
    main(args)
