#!/usr/bin/env python3
"""Convert GenBank format to FASTA format."""

import argparse
from pathlib import Path
from Bio import SeqIO


def convert(input_file, output_file):
    # Convert GenBank to FASTA
    records = list(SeqIO.parse(input_file, "genbank"))

    if not records:
        print(f"ERROR: No sequences found in {input_file}")
        exit(1)

    # Write sequences to FASTA format
    with open(output_file, 'w') as f:
        SeqIO.write(records, f, "fasta")

    print(f"Converted {input_file} -> {output_file}")
    print(f"  Sequences: {len(records)}")
    total_bp = sum(len(r.seq) for r in records)
    print(f"  Total length: {total_bp:,} bp")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GenBank to FASTA format")
    parser.add_argument("input_file", help="Input GenBank file")
    parser.add_argument("output_file", help="Output FASTA file")
    args = parser.parse_args()

    convert(args.input_file, args.output_file)
