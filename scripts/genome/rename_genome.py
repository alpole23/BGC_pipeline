#!/usr/bin/env python3
"""Rename a genome file based on the name map."""

import argparse
import json
import shutil
import sys
from pathlib import Path


def rename_genome(assembly_id, genome_file, name_map_file, output_dir="."):
    with open(name_map_file, 'r') as f:
        name_map = json.load(f)

    genome_path = Path(genome_file)

    if assembly_id not in name_map:
        print(f"ERROR: Assembly ID {assembly_id} not found in name map")
        sys.exit(1)

    new_name = name_map[assembly_id].replace(' ', '_')
    output_file = Path(output_dir) / f"{new_name}.gbff"

    shutil.copy2(genome_path, output_file)

    if not output_file.exists():
        print(f"ERROR: Output file was not created")
        sys.exit(1)

    file_size = output_file.stat().st_size
    if file_size == 0:
        print(f"ERROR: Output file is empty")
        sys.exit(1)

    print(f"Successfully created {output_file} ({file_size:,} bytes)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename genome file based on name map")
    parser.add_argument("assembly_id", help="Assembly ID (e.g., GCF_000001234.1)")
    parser.add_argument("genome_file", help="Path to genome file")
    parser.add_argument("name_map", help="Path to name_map.json")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    args = parser.parse_args()

    rename_genome(args.assembly_id, args.genome_file, args.name_map, args.output_dir)
