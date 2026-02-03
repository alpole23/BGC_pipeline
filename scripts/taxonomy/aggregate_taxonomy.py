#!/usr/bin/env python3
"""Aggregate BGC statistics by taxonomy into a hierarchical tree."""

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def build_tree(bgc_counts, taxonomy_map, genome_to_assembly):
    """Build hierarchical tree from taxonomy map and BGC counts."""
    tree = {
        'children': {},
        'genomes': [],
        'stats': defaultdict(lambda: defaultdict(int))
    }

    # Process each genome
    for genome_name, counts in bgc_counts.items():
        assembly_id = genome_to_assembly.get(genome_name)

        if not assembly_id or assembly_id not in taxonomy_map:
            print(f"Warning: Could not find taxonomy for genome {genome_name}")
            continue

        tax_info = taxonomy_map[assembly_id]
        lineage = tax_info['lineage']

        # Build path through tree
        current = tree
        path = []

        # Traverse taxonomy levels in order
        for rank in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            if rank not in lineage:
                continue

            tax_node = lineage[rank]
            tax_name = tax_node.get('name', '')
            tax_id = tax_node.get('id', '')

            if not tax_name:
                continue

            path.append(tax_name)

            # Create node if it doesn't exist
            if tax_name not in current['children']:
                current['children'][tax_name] = {
                    'name': tax_name,
                    'taxId': tax_id,
                    'rank': rank,
                    'parent_path': path[:-1].copy(),
                    'children': {},
                    'genomes': [],
                    'stats': defaultdict(lambda: defaultdict(int))
                }

            current = current['children'][tax_name]

        # Add genome to leaf node
        total_bgcs = sum(counts.values())
        genome_entry = {
            'assembly_id': assembly_id,
            'name': genome_name,
            'strain': tax_info.get('strain', ''),
            'total_bgcs': total_bgcs,
            'bgc_types': counts
        }
        current['genomes'].append(genome_entry)

    return tree


def aggregate_stats(node):
    """Recursively aggregate statistics from children and genomes."""

    # First, recursively process all children
    for child_name, child_node in node['children'].items():
        aggregate_stats(child_node)

    # Initialize aggregation variables
    all_genome_bgc_counts = []  # List of total BGC counts per genome
    bgc_type_stats = defaultdict(lambda: {'genome_count': 0, 'total_count': 0})
    genome_count = 0
    genomes_with_bgcs = 0
    genomes_without_bgcs = 0

    # Aggregate from direct genomes at this node
    for genome in node['genomes']:
        genome_count += 1
        total_bgcs = genome['total_bgcs']
        all_genome_bgc_counts.append(total_bgcs)

        if total_bgcs > 0:
            genomes_with_bgcs += 1
        else:
            genomes_without_bgcs += 1

        # Aggregate BGC type counts
        for bgc_type, count in genome['bgc_types'].items():
            if count > 0:
                bgc_type_stats[bgc_type]['genome_count'] += 1
                bgc_type_stats[bgc_type]['total_count'] += count

    # Aggregate from children
    for child_name, child_node in node['children'].items():
        child_stats = child_node['stats']

        genome_count += child_stats['genome_count']
        genomes_with_bgcs += child_stats['genomes_with_bgcs']
        genomes_without_bgcs += child_stats['genomes_without_bgcs']

        # Merge BGC type stats
        for bgc_type, type_stats in child_stats.get('bgc_type_distribution', {}).items():
            bgc_type_stats[bgc_type]['genome_count'] += type_stats['genome_count']
            bgc_type_stats[bgc_type]['total_count'] += type_stats['total_count']

        # Collect genome BGC counts from children
        if 'genome_bgc_counts' in child_stats:
            all_genome_bgc_counts.extend(child_stats['genome_bgc_counts'])

    # Calculate summary statistics
    total_bgcs = sum(all_genome_bgc_counts)
    avg_bgcs = total_bgcs / genome_count if genome_count > 0 else 0

    # Calculate standard deviation
    if genome_count > 1:
        variance = sum((x - avg_bgcs) ** 2 for x in all_genome_bgc_counts) / (genome_count - 1)
        std_bgcs = math.sqrt(variance)
    else:
        std_bgcs = 0

    # Build BGC type distribution with percentages
    bgc_distribution = {}
    for bgc_type, type_stats in bgc_type_stats.items():
        bgc_distribution[bgc_type] = {
            'total_count': type_stats['total_count'],
            'genome_count': type_stats['genome_count'],
            'avg_per_genome': type_stats['total_count'] / genome_count if genome_count > 0 else 0,
            'percentage_of_genomes': (type_stats['genome_count'] / genome_count * 100) if genome_count > 0 else 0
        }

    # Sort by total count descending
    bgc_distribution = dict(sorted(bgc_distribution.items(),
                                  key=lambda x: x[1]['total_count'],
                                  reverse=True))

    # Store aggregated stats
    node['stats'] = {
        'genome_count': genome_count,
        'total_bgcs': total_bgcs,
        'avg_bgcs_per_genome': round(avg_bgcs, 2),
        'std_bgcs_per_genome': round(std_bgcs, 2),
        'genomes_with_bgcs': genomes_with_bgcs,
        'genomes_without_bgcs': genomes_without_bgcs,
        'bgc_type_distribution': bgc_distribution,
        'genome_bgc_counts': all_genome_bgc_counts  # For parent aggregation
    }

    # Sort children by genome count (most diverse first)
    if node['children']:
        node['children'] = dict(sorted(node['children'].items(),
                                      key=lambda x: x[1]['stats']['genome_count'],
                                      reverse=True))


def cleanup(node):
    """Remove internal data used for aggregation."""
    if 'genome_bgc_counts' in node['stats']:
        del node['stats']['genome_bgc_counts']
    for child in node['children'].values():
        cleanup(child)


def main(taxonomy_map_file, region_counts_file, name_map_file, output_file, taxon):
    print("Loading taxonomy map...")
    with open(taxonomy_map_file, 'r') as f:
        taxonomy_map = json.load(f)

    print(f"Loaded {len(taxonomy_map)} assemblies")

    print("Loading BGC counts...")
    # Read region counts - maps genome name to BGC type counts
    bgc_counts = {}
    with open(region_counts_file, 'r') as f:
        # Skip comment lines
        lines = [line for line in f if not line.startswith('#')]
        reader = csv.DictReader(lines, delimiter='\t')
        for row in reader:
            genome_name = row['record']
            # Remove file extension
            if genome_name.endswith('.gbff'):
                genome_name = genome_name[:-5]
            # Remove non-count columns
            counts = {k: int(v) for k, v in row.items()
                     if k not in ['record', 'total_count', 'description']}
            bgc_counts[genome_name] = counts

    print(f"Loaded BGC counts for {len(bgc_counts)} genomes")

    # Load name map (assembly_id -> renamed_genome_name)
    print("Loading name map...")
    with open(name_map_file, 'r') as f:
        name_map = json.load(f)

    print(f"Loaded {len(name_map)} assembly-to-name mappings")

    # Create genome_to_assembly map (reversed name_map)
    genome_to_assembly = {}
    for assembly_id, genome_name in name_map.items():
        if genome_name in bgc_counts:
            genome_to_assembly[genome_name] = assembly_id

    print(f"Matched {len(genome_to_assembly)} genomes to assemblies")

    if len(genome_to_assembly) < len(bgc_counts):
        unmatched = set(bgc_counts.keys()) - set(genome_to_assembly.keys())
        print(f"Warning: {len(unmatched)} genomes could not be matched to assemblies")
        if len(unmatched) <= 10:
            for genome_name in list(unmatched)[:10]:
                print(f"  Unmatched: {genome_name}")

    # Build hierarchical tree structure
    print("Building taxonomic tree...")
    tree = build_tree(bgc_counts, taxonomy_map, genome_to_assembly)

    # Aggregate statistics bottom-up
    print("Aggregating statistics...")
    aggregate_stats(tree)

    # Clean up internal data
    cleanup(tree)

    # Find the root node (should be at domain level)
    root_node = None
    if tree['children']:
        # Get the first (and typically only) domain
        root_name = next(iter(tree['children']))
        root_node = tree['children'][root_name]
    else:
        # No children, create a dummy root
        root_node = {
            'name': taxon,
            'taxId': 0,
            'rank': 'unknown',
            'stats': tree['stats'],
            'children': {},
            'genomes': tree['genomes']
        }

    # Build output structure
    output = {
        'metadata': {
            'pipeline_version': '1.0',
            'query_taxon': taxon,
            'total_genomes': root_node['stats']['genome_count'],
            'taxon_levels': ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        },
        'tree': root_node
    }

    # Write output
    print("Writing taxonomy_tree.json...")
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    print("Done!")
    print(f"Total genomes in tree: {output['metadata']['total_genomes']}")
    print(f"Total BGCs: {root_node['stats']['total_bgcs']}")
    print(f"Root node: {root_node['name']} ({root_node['rank']})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate BGC statistics by taxonomy")
    parser.add_argument("taxonomy_map", help="Path to taxonomy_map.json")
    parser.add_argument("region_counts", help="Path to region_counts.tsv")
    parser.add_argument("name_map", help="Path to name_map.json")
    parser.add_argument("output_file", help="Path to output JSON file")
    parser.add_argument("--taxon", default="Unknown", help="Taxon name for metadata")
    args = parser.parse_args()

    main(args.taxonomy_map, args.region_counts, args.name_map, args.output_file, args.taxon)
