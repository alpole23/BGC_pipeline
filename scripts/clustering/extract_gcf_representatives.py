#!/usr/bin/env python3
"""Extract GCF representative data from BiG-SCAPE database and antiSMASH results."""

import argparse
import json
import os
import sqlite3
import sys
from pathlib import Path

import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.constants import GENE_COLORS, GCF_COLORS
from utils.antismash_parser import (
    parse_location, get_gene_function_category, get_gene_color,
    find_antismash_json, extract_domain_hits, load_antismash_json,
    extract_genes_from_region, find_record_index, build_record_index_map
)
from utils.gene_diagram import generate_gene_svg


def load_kcb_lookup(tabulation_file):
    """Load KCB hit data from tabulation file into a lookup dict.

    Returns: dict mapping (genome, region_name) -> {'kcb_hit': str, 'kcb_acc': str}

    Note: Uses region_name (e.g., "40.1") as the unique identifier since it
    combines record_index and region number, avoiding collisions.
    """
    kcb_lookup = {}
    if not tabulation_file or not os.path.exists(tabulation_file):
        return kcb_lookup

    try:
        df = pd.read_csv(tabulation_file, sep='\t', comment='#')
        for _, row in df.iterrows():
            genome = row.get('file', '')
            region_name = row.get('region_name', '')
            if genome and region_name:
                # Use (genome, region_name) as key - region_name is unique within genome
                key = (genome, str(region_name))
                # Handle NaN values from pandas - convert to empty string
                kcb_hit = row.get('KCB_hit', '')
                kcb_acc = row.get('KCB_acc', '')
                if pd.isna(kcb_hit):
                    kcb_hit = ''
                if pd.isna(kcb_acc):
                    kcb_acc = ''
                kcb_lookup[key] = {
                    'kcb_hit': str(kcb_hit) if kcb_hit else '',
                    'kcb_acc': str(kcb_acc) if kcb_acc else ''
                }
    except Exception as e:
        print(f"Warning: Could not load KCB data from {tabulation_file}: {e}")

    return kcb_lookup


def extract_genome_gcf_mapping(cursor, antismash_dir):
    """Extract mapping of genome names to their GCF families.

    Returns:
        genome_to_gcfs: {genome_name: [family_id, ...]}
        family_metadata: {family_id: {product, member_count, color}}
        bgc_to_gcf: {(genome_name, region_name): {family_id, member_count}}

    Note: bgc_to_gcf uses region_name (e.g., "40.1") as the unique identifier.
    """
    genome_to_gcfs = {}
    family_metadata = {}
    bgc_to_gcf = {}  # Maps (genome, region_name) -> {family_id, member_count}

    # Build record_id -> record_index mapping from antiSMASH JSONs
    record_index_map = build_record_index_map(antismash_dir)

    # Get all families with member counts
    cursor.execute('''
        SELECT
            f.id as family_id,
            f.cutoff,
            f.bin_label,
            br.product,
            COUNT(brf.record_id) as member_count
        FROM family f
        JOIN bgc_record br ON f.center_id = br.id
        LEFT JOIN bgc_record_family brf ON f.id = brf.family_id
        GROUP BY f.id
        ORDER BY member_count DESC
    ''')

    families = cursor.fetchall()
    for i, family in enumerate(families):
        family_id = family['family_id']
        color = GCF_COLORS[i % len(GCF_COLORS)]
        family_metadata[family_id] = {
            'product': family['product'] or 'unknown',
            'member_count': family['member_count'],
            'color': color
        }

    # Get all BGC-to-family mappings with genome paths and region numbers
    cursor.execute('''
        SELECT
            brf.family_id,
            g.path as gbk_path,
            br.record_number as region_number
        FROM bgc_record_family brf
        JOIN bgc_record br ON brf.record_id = br.id
        JOIN gbk g ON br.gbk_id = g.id
    ''')

    rows = cursor.fetchall()
    for row in rows:
        gbk_path = row['gbk_path']
        family_id = row['family_id']
        region_number = row['region_number']

        # Extract genome name from path
        # BiG-SCAPE stores paths like: .../antismash_input/Pantoea_ananatis_19-20/...
        path_parts = Path(gbk_path).parts
        genome_name = None
        for idx, part in enumerate(path_parts):
            if part == 'antismash_input' and idx + 1 < len(path_parts):
                genome_name = path_parts[idx + 1]
                break

        if not genome_name:
            # Fallback: try to find from directory structure
            for part in reversed(path_parts[:-1]):  # Exclude filename
                if part and not part.startswith('.') and 'region' not in part:
                    genome_name = part
                    break

        # Extract record_id from GBK filename (e.g., "VNIC01000040.1.region001.gbk" -> "VNIC01000040")
        gbk_basename = Path(gbk_path).stem
        record_id = gbk_basename.split('.region')[0] if '.region' in gbk_basename else gbk_basename
        # Remove trailing .1 if present (NCBI version suffix)
        if record_id.endswith('.1'):
            record_id = record_id[:-2]

        if genome_name:
            if genome_name not in genome_to_gcfs:
                genome_to_gcfs[genome_name] = set()
            genome_to_gcfs[genome_name].add(family_id)

            # Build BGC-to-GCF mapping using region_name for unique identification
            if region_number and record_id:
                # Look up record_index to construct region_name
                record_index = record_index_map.get((genome_name, record_id), 0)
                if record_index > 0:
                    region_name = f"{record_index}.{region_number}"
                    meta = family_metadata.get(family_id, {})
                    bgc_to_gcf[(genome_name, region_name)] = {
                        'family_id': family_id,
                        'member_count': meta.get('member_count', 1)
                    }

    # Convert sets to sorted lists
    genome_to_gcfs = {k: sorted(v) for k, v in genome_to_gcfs.items()}

    return genome_to_gcfs, family_metadata, bgc_to_gcf


def extract_gcf_representatives(bigscape_dir, antismash_dir, output_file, taxon=None, tabulation_file=None):
    """Extract GCF representative data from BiG-SCAPE database."""

    # Use provided taxon or fall back to extracting from antismash_dir
    taxon_name = taxon if taxon else Path(antismash_dir).name

    # Load KCB lookup for novelty detection
    kcb_lookup = load_kcb_lookup(tabulation_file)
    if kcb_lookup:
        print(f"Loaded KCB data for {len(kcb_lookup)} BGC regions")

    # Find SQLite database
    db_files = list(Path(bigscape_dir).glob("*.db"))
    if not db_files:
        print("Error: No BiG-SCAPE database found")
        result = {"error": "No BiG-SCAPE database found", "gcfs": [], "summary": {},
                  "genome_gcf_mapping": {}, "family_metadata": {}, "bgc_to_gcf": {}}
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        return

    db_path = db_files[0]
    print(f"Using database: {db_path}")
    print(f"Taxon: {taxon_name}")

    gcfs = []
    genome_gcf_mapping = {}
    family_metadata = {}
    bgc_to_gcf = {}

    try:
        conn = sqlite3.connect(str(db_path))
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        # Extract genome-to-GCF mapping for tree coloring and BGC-to-GCF for novel BGCs
        genome_gcf_mapping, family_metadata, bgc_to_gcf = extract_genome_gcf_mapping(cursor, antismash_dir)
        print(f"Mapped {len(genome_gcf_mapping)} genomes to GCFs")
        print(f"Mapped {len(bgc_to_gcf)} individual BGCs to GCF families")

        # Get all families with their representatives and member counts
        cursor.execute('''
            SELECT
                f.id as family_id,
                f.center_id,
                f.cutoff,
                f.bin_label,
                COUNT(brf.record_id) as member_count
            FROM family f
            LEFT JOIN bgc_record_family brf ON f.id = brf.family_id
            GROUP BY f.id
            ORDER BY member_count DESC
        ''')

        families = cursor.fetchall()
        print(f"Found {len(families)} gene cluster families")

        for family in families:
            family_id = family['family_id']
            center_id = family['center_id']
            member_count = family['member_count']
            is_singleton = member_count == 1

            # Get representative BGC details
            cursor.execute('''
                SELECT
                    br.id,
                    br.record_number,
                    g.path,
                    g.organism,
                    br.product,
                    br.category
                FROM bgc_record br
                JOIN gbk g ON br.gbk_id = g.id
                WHERE br.id = ?
            ''', (center_id,))

            rep_row = cursor.fetchone()
            if not rep_row:
                print(f"Warning: Could not find representative for family {family_id}")
                continue

            gbk_path = rep_row['path']
            organism = rep_row['organism'] or 'Unknown organism'
            product = rep_row['product'] or 'Unknown product'
            category = rep_row['category'] or 'other'
            region_number = rep_row['record_number']

            # Extract genome name from path
            # BiG-SCAPE stores paths like: .../antismash_input/Pantoea_ananatis_19-20/...
            path_parts = Path(gbk_path).parts
            genome_name = None
            for part in path_parts:
                if 'antismash_input' in path_parts:
                    idx = path_parts.index('antismash_input')
                    if idx + 1 < len(path_parts):
                        genome_name = path_parts[idx + 1]
                        break

            if not genome_name:
                # Try to extract from the GBK filename
                gbk_name = Path(gbk_path).stem
                # Region GBK files are like: JABDZE010000001.1.region001
                # Main dir might be genome name
                genome_name = gbk_name.split('.region')[0] if '.region' in gbk_name else gbk_name

            # Also extract record_id from GBK filename
            gbk_basename = Path(gbk_path).stem
            record_id = gbk_basename.split('.region')[0] if '.region' in gbk_basename else gbk_basename
            # Remove trailing .1 if present
            if record_id.endswith('.1'):
                record_id = record_id[:-2]

            # Find antiSMASH JSON and extract genes
            antismash_json_path = find_antismash_json(genome_name, antismash_dir)
            antismash_data = load_antismash_json(antismash_json_path) if antismash_json_path else None
            genes = []
            region_start = 0
            region_end = 0
            svg_diagram = ''

            record_index = 1  # Default
            enhanced_analysis = None
            if antismash_data:
                genes, region_start, region_end, parsed_product, enhanced_analysis = extract_genes_from_region(
                    antismash_data, record_id, region_number
                )
                if parsed_product:
                    product = parsed_product

                if genes:
                    svg_diagram = generate_gene_svg(genes, region_start, region_end)

                # Find the correct record index for the antiSMASH anchor
                record_index = find_record_index(antismash_data, record_id)
            else:
                print(f"Warning: Could not find antiSMASH JSON for genome {genome_name}")

            # Build antiSMASH link
            # Link is relative from main_data_visualization folder to antismash_results
            # Path: main_analysis_results/{taxon}/main_data_visualization/ -> antismash_results/{taxon}/{genome}/
            # Anchor format: #r{record_index}c{region_number}
            antismash_link = f"../../../antismash_results/{taxon_name}/{genome_name}/index.html#r{record_index}c{region_number}" if genome_name else None

            # Look up KCB hit for this representative BGC
            # Key format: (genome, region_name) where region_name = "{record_index}.{region_number}"
            region_name = f"{record_index}.{region_number}"
            kcb_info = kcb_lookup.get((genome_name, region_name), {})
            kcb_hit = kcb_info.get('kcb_hit', '')
            kcb_acc = kcb_info.get('kcb_acc', '')

            gcf_data = {
                'family_id': family_id,
                'member_count': member_count,
                'is_singleton': is_singleton,
                'cutoff': family['cutoff'],
                'product': product,
                'category': category,
                'organism': organism,
                'genome_name': genome_name or 'unknown',
                'record_id': record_id,
                'region_number': region_number,
                'region_start': region_start,
                'region_end': region_end,
                'genes': genes,
                'svg_diagram': svg_diagram,
                'antismash_link': antismash_link,
                'enhanced_analysis': enhanced_analysis,
                'kcb_hit': kcb_hit,
                'kcb_acc': kcb_acc
            }

            gcfs.append(gcf_data)

        conn.close()

    except sqlite3.Error as e:
        print(f"Database error: {e}")
        result = {"error": str(e), "gcfs": [], "summary": {},
                  "genome_gcf_mapping": {}, "family_metadata": {}, "bgc_to_gcf": {}}
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        return

    # Calculate summary
    total = len(gcfs)
    singletons = sum(1 for g in gcfs if g['is_singleton'])
    clusters = total - singletons

    # Convert bgc_to_gcf keys from tuples to strings for JSON serialization
    # Key format: "genome|region_name" (e.g., "genome|40.1")
    bgc_to_gcf_serializable = {
        f"{genome}|{region_name}": info
        for (genome, region_name), info in bgc_to_gcf.items()
    }

    result = {
        'gcfs': gcfs,
        'summary': {
            'total': total,
            'singletons': singletons,
            'clusters': clusters
        },
        'genome_gcf_mapping': genome_gcf_mapping,
        'family_metadata': family_metadata,
        'bgc_to_gcf': bgc_to_gcf_serializable
    }

    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"Extracted {total} GCFs ({clusters} clusters, {singletons} singletons)")
    print(f"Genome-to-GCF mapping: {len(genome_gcf_mapping)} genomes")
    print(f"Output written to: {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract GCF representative data from BiG-SCAPE')
    parser.add_argument('bigscape_dir', help='Path to BiG-SCAPE output directory')
    parser.add_argument('antismash_dir', help='Path to antiSMASH results directory')
    parser.add_argument('output_file', help='Path to output JSON file')
    parser.add_argument('--taxon', help='Taxon name for antiSMASH link paths (optional)')
    parser.add_argument('--tabulation', help='Path to region_tabulation.tsv for KCB hit lookup (optional)')

    args = parser.parse_args()

    extract_gcf_representatives(args.bigscape_dir, args.antismash_dir, args.output_file, args.taxon, args.tabulation)
