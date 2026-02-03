#!/usr/bin/env python3
"""Create a name map from assembly info to clean genome names."""

import argparse
import json
from collections import Counter
from pathlib import Path


CHAR_MAP = {
    ' ': '_',
    '\'': '_',
    '/': '_',
    '[': '',
    ']': '',
    '(': '',
    ')': '',
    ',': '_',
    '=': '_'
}


def clean_name(name):
    """Clean a name by replacing special characters."""
    return ''.join(CHAR_MAP.get(c, c) for c in name)


def create_name_map(assembly_info_file, output_file):
    name_map = {}
    # Track discovered names case-insensitively to avoid GTDB-Tk collisions
    # GTDB-Tk normalizes genome names, so MDCuke and MDcuke would collide
    discovered_names_lower = Counter()

    with open(assembly_info_file, 'r') as f:
        for line in f:
            if line.startswith('GCF') or line.startswith('GCA'):
                fields = line.strip().split('\t')

                # Extract fields
                assembly_id = fields[0]
                organism_name = fields[1] if len(fields) > 1 else ''
                strain = fields[2] if len(fields) > 2 else ''
                isolate = fields[4] if len(fields) > 4 else ''

                # Clean organism name
                org_name_clean = clean_name(organism_name)

                # Clean strain name
                strain_clean = clean_name(strain)

                # Clean isolate name
                isolate_clean = clean_name(isolate)

                # Build the base name with fallback logic:
                # 1. Use strain if available
                # 2. Use isolate if strain is empty
                # 3. Use assembly_id if both are empty
                if strain_clean and org_name_clean.endswith(strain_clean):
                    # Strain is redundant (already in organism name)
                    base_name = org_name_clean
                elif strain_clean:
                    # Strain is different, append it
                    base_name = f"{org_name_clean}_{strain_clean}"
                elif isolate_clean:
                    # No strain but isolate available, use isolate
                    base_name = f"{org_name_clean}_{isolate_clean}"
                else:
                    # No strain or isolate info - use assembly ID
                    base_name = f"{org_name_clean}_{assembly_id}"

                # Handle duplicates case-insensitively (GTDB-Tk normalizes names)
                base_name_lower = base_name.lower()
                discovered_names_lower[base_name_lower] += 1
                occurrences = discovered_names_lower[base_name_lower]

                if occurrences == 1:
                    final_name = base_name
                else:
                    # Add counter for duplicates (e.g., _2, _3, etc.)
                    final_name = f"{base_name}_{occurrences}"

                name_map[assembly_id] = final_name

    with open(output_file, 'w') as f:
        json.dump(name_map, f, indent=2)

    # Print summary
    print(f"Created name map for {len(name_map)} assemblies")
    print(f"Unique base names (case-insensitive): {len(discovered_names_lower)}")
    print(f"Names with duplicates: {sum(1 for count in discovered_names_lower.values() if count > 1)}")
    print("\nFirst 10 mappings:")
    for i, (assembly_id, name) in enumerate(list(name_map.items())[:10]):
        print(f"  {assembly_id} -> {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create name map from assembly info")
    parser.add_argument("assembly_info", help="Path to assembly info TSV file")
    parser.add_argument("output_file", help="Path to output JSON file")
    args = parser.parse_args()

    create_name_map(args.assembly_info, args.output_file)
