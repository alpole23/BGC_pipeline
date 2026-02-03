#!/usr/bin/env python3
"""Extract statistics from BiG-SCAPE clustering results."""

import argparse
import json
import os
from collections import defaultdict


def extract_stats(bigscape_dir, output_file):
    stats = {}

    # Find the output directory (contains timestamp)
    output_dir = os.path.join(bigscape_dir, "output_files")
    if not os.path.exists(output_dir):
        stats["error"] = "BiG-SCAPE output directory not found"
        with open(output_file, "w") as f:
            json.dump(stats, f, indent=2)
        return

    # Find all timestamped directories
    timestamp_dirs = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]

    if not timestamp_dirs:
        stats["error"] = "No BiG-SCAPE results directories found"
        with open(output_file, "w") as f:
            json.dump(stats, f, indent=2)
        return

    # Process each cutoff (there may be multiple)
    cutoff_stats = {}

    for timestamp_dir in timestamp_dirs:
        # Extract cutoff value from directory name (e.g., 2026-01-04_19-41-04_c0.5 -> 0.5)
        cutoff_str = timestamp_dir.split('_c')[-1] if '_c' in timestamp_dir else 'unknown'

        timestamp_path = os.path.join(output_dir, timestamp_dir)

        # Find all BGC class directories
        bgc_classes = [d for d in os.listdir(timestamp_path)
                      if os.path.isdir(os.path.join(timestamp_path, d))
                      and d != 'record_annotations']

        # Collect data across all BGC classes
        all_families = set()
        all_bgcs = set()
        family_sizes = defaultdict(int)
        class_families = defaultdict(set)
        class_bgcs = defaultdict(set)
        mibig_families = set()
        mibig_bgcs_found = False  # Track if ANY MIBiG BGCs exist in the dataset

        for bgc_class in bgc_classes:
            class_dir = os.path.join(timestamp_path, bgc_class)
            clustering_file = os.path.join(class_dir, f"{bgc_class}_clustering_c{cutoff_str}.tsv")

            if not os.path.exists(clustering_file):
                continue

            # Parse clustering file
            with open(clustering_file, 'r') as f:
                header = f.readline().strip().split('\t')

                # Find column indices
                try:
                    family_idx = header.index('Family')
                    record_idx = header.index('Record')
                except ValueError:
                    continue

                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) <= max(family_idx, record_idx):
                        continue

                    family = parts[family_idx]
                    record = parts[record_idx]

                    # Track overall families and BGCs
                    all_families.add(family)
                    all_bgcs.add(record)
                    family_sizes[family] += 1

                    # Track per-class statistics
                    class_families[bgc_class].add(family)
                    class_bgcs[bgc_class].add(record)

                    # Check for MIBiG references (MIBiG accessions start with 'BGC')
                    if 'BGC' in record or 'mibig' in record.lower():
                        mibig_families.add(family)
                        mibig_bgcs_found = True

        # Calculate statistics for this cutoff
        if all_families:
            family_size_list = list(family_sizes.values())

            cutoff_data = {
                "cutoff": float(cutoff_str) if cutoff_str != 'unknown' else cutoff_str,
                "total_families": len(all_families),
                "total_bgcs": len(all_bgcs),
                "avg_bgcs_per_family": round(sum(family_size_list) / len(family_size_list), 2) if family_size_list else 0,
                "max_bgcs_per_family": max(family_size_list) if family_size_list else 0,
                "min_bgcs_per_family": min(family_size_list) if family_size_list else 0,
                "singleton_families": sum(1 for size in family_size_list if size == 1),
                "mibig_included": mibig_bgcs_found,
                "families_with_mibig": len(mibig_families) if mibig_bgcs_found else None,
                "bgc_classes": {}
            }

            # Add per-class breakdown
            for bgc_class in sorted(class_families.keys()):
                cutoff_data["bgc_classes"][bgc_class] = {
                    "families": len(class_families[bgc_class]),
                    "bgcs": len(class_bgcs[bgc_class])
                }

            cutoff_stats[cutoff_str] = cutoff_data

    # Store all cutoff statistics
    stats["cutoffs"] = cutoff_stats

    # Add highest cutoff to top level for easy access (when multiple cutoffs exist)
    if cutoff_stats:
        # Sort cutoffs by value (highest first)
        sorted_cutoffs = sorted(cutoff_stats.items(),
                                key=lambda x: float(x[0]) if x[0] != 'unknown' else 0,
                                reverse=True)
        highest_cutoff = sorted_cutoffs[0][1]
        stats.update(highest_cutoff)

    # Write statistics to JSON file
    with open(output_file, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"BiG-SCAPE statistics extracted successfully")
    if stats.get("total_families"):
        print(f"Total Families: {stats['total_families']}")
        print(f"Total BGCs: {stats['total_bgcs']}")
        print(f"Avg BGCs per Family: {stats['avg_bgcs_per_family']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract BiG-SCAPE clustering statistics")
    parser.add_argument("bigscape_dir", help="Path to BiG-SCAPE output directory")
    parser.add_argument("output_file", help="Path to output JSON file")
    args = parser.parse_args()

    extract_stats(args.bigscape_dir, args.output_file)
