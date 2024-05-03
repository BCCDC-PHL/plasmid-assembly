#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys

def parse_ale(ale_path):
    """
    """
    ale = {}
    field_name_translation = {
        'ALE_score': 'ale_score',
        'numContigs': 'num_contigs',
        'totalAssemLen': 'total_assembly_length',
        'placeAvg': 'place_avg',
        'insertAvg': 'insert_avg',
        'kmerAvg': 'kmer_avg',
        'depthScoreAvg': 'depth_score_avg',
        'depthAvg': 'depth_avg',
        'totalReads': 'total_reads',
        'totalMappedReads': 'total_mapped_reads',
        'totalUnmappedReads': 'total_unmapped_reads',
        'totalPlacedReads': 'total_placed_reads',
        'readAvgLen': 'read_avg_length',
        'avgReadOverlap': 'average_read_overlap',
        'Reference': 'reference',
    }
    float_fields = [
        'ale_score',
        'place_avg',
        'insert_avg',
        'kmer_avg',
        'depth_score_avg',
        'depth_avg',
        'total_mapped_reads',
        'total_unmapped_reads',
        'total_placed_reads',
        'read_avg_length',
        'average_read_overlap',
    ]
    int_fields = [
        'num_contigs',
        'total_assembly_length',
        'total_reads',
        'total_mapped_reads',
        'total_unmapped_reads',
        'total_placed_reads',
    ]
    with open(ale_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") and not line.startswith("# contig"):
                line = line.strip("#").strip()
                line_split = line.split(": ")
                field_name = line_split[0]
                value = line_split[1]
                translated_field_name = field_name_translation.get(field_name, field_name)
                if translated_field_name in float_fields:
                    try:
                        value = float(value)
                    except ValueError as e:
                        value = None
                elif translated_field_name in int_fields:
                    try:
                        value = int(value)
                    except ValueError as e:
                        value = None
                ale[translated_field_name] = value
            else:
                break

    return ale

    
def main(args):
    output_fieldnames = [
        'ale_score',
        'num_contigs',
        'total_assembly_length',
        'place_avg',
        'insert_avg',
        'kmer_avg',
        'depth_score_avg',
        'depth_avg',
        'total_reads',
        'total_mapped_reads',
        'total_unmapped_reads',
        'total_placed_reads',
        'read_avg_length',
        'average_read_overlap',
    ]
    ale = parse_ale(args.input)
    if args.sample_id:
        ale["sample_id"] = args.sample_id
        output_fieldnames.insert(0, 'sample_id')

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='unix', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()
    writer.writerow(ale)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input file")
    parser.add_argument('-s', '--sample-id', help='sample ID')
    args = parser.parse_args()
    main(args)
