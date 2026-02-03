#!/usr/bin/env python3
"""
BGC Analysis Visualization Script

Generates interactive HTML reports and visualizations for biosynthetic gene cluster analysis.
Includes:
- Tab-based HTML report with sections (Overview, Taxonomy, BGC Distribution, Genomes, Clustering)
- Summary statistics dashboard with tabulation stats
- Taxonomic distribution tree (interactive expandable)
- GCF x Taxonomy heatmap and distribution analysis (replaces tree visualization)
- Donut chart for BGC type distribution
- BGC histogram showing distribution across genomes
- Searchable genome table
- Individual genome metadata pages
- Clustering statistics (BiG-SCAPE and BiG-SLiCE)
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Wedge, FancyBboxPatch
import seaborn as sns
from pathlib import Path
import json
import numpy as np
import shutil
import re
import sys
import urllib.parse
from collections import defaultdict

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from utils.parsers import (
    parse_trace_file, parse_duration, parse_memory,
    format_bytes, format_duration_str
)
from utils.trace import (
    aggregate_trace_by_process, generate_resource_usage_html
)

# =============================================================================
# BGC CLASS COLOR PALETTE
# =============================================================================
# Consistent color scheme inspired by antiSMASH conventions
# These colors are used across all visualizations for easy pattern recognition

BGC_COLORS = {
    # PKS types - Blue family
    'T1PKS': '#1f77b4',      # Strong blue
    'T2PKS': '#4a90d9',      # Medium blue
    'T3PKS': '#7eb3ed',      # Light blue
    'transAT-PKS': '#2c5aa0', # Dark blue
    'PKS-like': '#5dade2',   # Sky blue
    'hglE-KS': '#85c1e9',    # Pale blue

    # NRPS types - Red/Orange family
    'NRPS': '#d62728',       # Strong red
    'NRPS-like': '#ff6b6b',  # Light red
    'thioamide-NRP': '#e74c3c', # Tomato red
    'NAPAA': '#c0392b',      # Dark red

    # Hybrid types - Purple family
    'PKS-NRPS_Hybrids': '#9b59b6',  # Purple (for aggregated hybrids)
    'NRPS-PKS_Hybrids': '#8e44ad',  # Dark purple

    # RiPPs - Green family
    'RiPP': '#27ae60',       # Strong green
    'lanthipeptide': '#2ecc71',     # Emerald
    'lanthipeptide-class-i': '#27ae60',
    'lanthipeptide-class-ii': '#229954',
    'lanthipeptide-class-iii': '#1e8449',
    'lanthipeptide-class-iv': '#196f3d',
    'lanthipeptide-class-v': '#145a32',
    'thiopeptide': '#58d68d',
    'LAP': '#82e0aa',
    'lassopeptide': '#abebc6',
    'sactipeptide': '#d5f5e3',
    'bottromycin': '#a9dfbf',
    'cyanobactin': '#73c6b6',
    'microviridin': '#45b39d',
    'proteusin': '#16a085',
    'RRE-containing': '#138d75',
    'fungal-RiPP': '#117a65',
    'ranthipeptide': '#0e6655',
    'redox-cofactor': '#0b5345',
    'RiPP-like': '#7dcea0',
    'thioamitides': '#52be80',
    'epipeptide': '#48c9b0',
    'guanidinotides': '#1abc9c',
    'glycocin': '#17a589',
    'triceptide': '#148f77',
    'spliceotide': '#117864',
    'methanobactin': '#0e6251',
    'cyclic-lactone-autoinducer': '#85929e',
    'darobactin': '#76d7c4',
    'rcdpeptide': '#45b39d',

    # Terpenes - Yellow/Gold family
    'terpene': '#f39c12',    # Orange-yellow

    # Saccharides - Brown family
    'saccharide': '#a0522d', # Sienna brown
    'oligosaccharide': '#cd853f',
    'polysaccharide': '#8b4513',
    'amglyccycl': '#d2691e',

    # Alkaloids/Other nitrogen - Teal family
    'alkaloid': '#008080',
    'indole': '#20b2aa',
    'NI-siderophore': '#40e0d0',

    # Siderophores - Cyan family
    'siderophore': '#00ced1',
    'NAGGN': '#00bfff',

    # Fatty acids / Lipids - Olive family
    'fatty_acid': '#808000',
    'PUFA': '#9acd32',
    'ladderane': '#6b8e23',
    'hserlactone': '#556b2f',
    'acyl_amino_acids': '#8fbc8f',
    'N-acyl amino acid': '#8fbc8f',

    # Phosphonates - Pink family
    'phosphonate': '#ff69b4',
    'phosphonate-like': '#ffb6c1',

    # Aromatic compounds - Magenta family
    'arylpolyene': '#ff00ff',
    'resorcinol': '#da70d6',
    'stilbene': '#ee82ee',
    'phenazine': '#dda0dd',
    'aminocoumarin': '#ba55d3',

    # Nucleosides - Gray family
    'nucleoside': '#708090',

    # Bacteriocins - Dark cyan
    'bacteriocin': '#008b8b',
    'RaS-RiPP': '#5f9ea0',

    # Beta-lactams - Coral
    'betalactam': '#ff7f50',
    'beta-lactam': '#ff7f50',

    # Ectoine - Light purple
    'ectoine': '#dda0dd',
    'ectoine-like': '#e6e6fa',

    # Melanin - Dark gray
    'melanin': '#2f4f4f',

    # Butyrolactone - Peach
    'butyrolactone': '#ffdab9',

    # Blactam - Salmon
    'blactam': '#fa8072',

    # CDPS - Lavender
    'CDPS': '#e6e6fa',

    # Furan - Wheat
    'furan': '#f5deb3',

    # Prodigiosin - Crimson
    'prodigiosin': '#dc143c',

    # Cyanide - Light steel blue
    'cyanide': '#b0c4de',
    'hydrogen-cyanide': '#b0c4de',

    # Linaridin - Medium purple
    'linaridin': '#9370db',
    'linear-azol(in)e-containing-peptide': '#9370db',

    # Opine - Thistle
    'opine-like-metallophore': '#d8bfd8',

    # Other/Unknown - Gray
    'other': '#95a5a6',
    'Other': '#95a5a6',
    'unknown': '#bdc3c7',
    'Unknown': '#bdc3c7',
    'NA': '#ecf0f1',
}

def get_bgc_color(bgc_type):
    '''Get color for a BGC type, handling hybrids and unknown types'''
    if bgc_type in BGC_COLORS:
        return BGC_COLORS[bgc_type]

    # Check if it's a hybrid (contains + or -)
    if '+' in bgc_type or '-' in bgc_type:
        # Try to match first component
        parts = bgc_type.replace('+', '-').split('-')
        for part in parts:
            if part in BGC_COLORS:
                return BGC_COLORS[part]
        return '#9b59b6'  # Default purple for hybrids

    # Check for partial matches
    bgc_lower = bgc_type.lower()
    for key, color in BGC_COLORS.items():
        if key.lower() in bgc_lower or bgc_lower in key.lower():
            return color

    return '#95a5a6'  # Default gray for unknown

def create_genome_metadata_pages(counts_file, assembly_info, name_map, outdir, taxon, taxonomy_map_data=None, tabulation_file=None):
    '''Create individual HTML pages for each genome with metadata and KCB data'''
    import re
    genome_dir = outdir / 'genomes'
    genome_dir.mkdir(exist_ok=True)

    # Clean taxon name for URL - match Nextflow sanitizeTaxon function
    taxon_clean = re.sub(r'[^a-zA-Z0-9_]', '_', taxon)
    taxon_clean = re.sub(r'_+', '_', taxon_clean).strip('_')

    # Read counts to get genome list
    counts_df = pd.read_csv(counts_file, sep='\t', comment='#')

    # Read assembly info
    assembly_df = pd.read_csv(assembly_info, sep='\t')

    # Load tabulation data for KCB information
    tab_df = None
    if tabulation_file:
        try:
            tab_df = pd.read_csv(tabulation_file, sep='\t')
            # Fill NaN values with empty strings for consistent filtering
            for col in ['KCB_hit', 'KCB_acc', 'KCB_sim']:
                if col in tab_df.columns:
                    tab_df[col] = tab_df[col].fillna('')
        except Exception as e:
            print(f"Warning: Could not load tabulation file: {e}")
            tab_df = None

    # Load name map and create reverse mapping (genome_name -> assembly_id)
    with open(name_map, 'r') as f:
        name_map_data = json.load(f)
    reverse_map = {v: k for k, v in name_map_data.items()}

    for _, row in counts_df.iterrows():
        genome_name = row['record'].replace('.gbff', '')

        # Look up assembly ID using reverse name map
        assembly_id = reverse_map.get(genome_name, None)

        # Find matching assembly metadata
        metadata_rows = assembly_df[assembly_df.iloc[:, 0] == assembly_id] if assembly_id else pd.DataFrame()

        if not metadata_rows.empty:
            metadata = metadata_rows.iloc[0]

            def get_value(idx):
                if len(metadata) > idx:
                    val = metadata.iloc[idx]
                    if pd.isna(val) or val == '' or str(val) == 'nan':
                        return 'N/A'
                    return str(val)
                return 'N/A'

            organism_name = get_value(1)
            strain = get_value(2)
            isolation_source = get_value(3)
            isolate = get_value(4)
            notes = get_value(5)
            checkm_completeness = get_value(6)
            checkm_contamination = get_value(7)
            organism_tax_id = get_value(8)

            # Display taxonomy lineage if available
            taxonomy_html = ''
            if taxonomy_map_data and assembly_id and assembly_id in taxonomy_map_data:
                genome_taxonomy = taxonomy_map_data[assembly_id]
                lineage = genome_taxonomy.get('lineage', {})

                lineage_parts = []
                for rank in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                    if rank in lineage:
                        tax_node = lineage[rank]
                        name = tax_node.get('name', 'N/A')
                        if name not in ['N/A', '', 'nan', 'None']:
                            lineage_parts.append(f"<strong>{rank.capitalize()}:</strong> {name}")

                if lineage_parts:
                    lineage_html = '<br>'.join(lineage_parts)
                    taxonomy_html = f'''
        <tr><td colspan="2"><strong>Taxonomic Lineage</strong></td></tr>
        <tr><td colspan="2">{lineage_html}</td></tr>'''

            assembly_table = f'''
    <table>
        <tr><th>Field</th><th>Value</th></tr>
        <tr><td>Assembly Accession</td><td>{assembly_id}</td></tr>
        <tr><td>Organism Tax ID</td><td>{organism_tax_id}</td></tr>
        <tr><td>Organism Name</td><td>{organism_name}</td></tr>
        <tr><td>Strain</td><td>{strain}</td></tr>
        <tr><td>Isolate</td><td>{isolate}</td></tr>
        <tr><td>Isolation Source</td><td>{isolation_source}</td></tr>
        <tr><td>Assembly Notes</td><td>{notes}</td></tr>
        <tr><td>CheckM Completeness</td><td>{checkm_completeness}</td></tr>
        <tr><td>CheckM Contamination</td><td>{checkm_contamination}</td></tr>
        {taxonomy_html}
    </table>'''
        else:
            assembly_table = '<p><em>Assembly metadata not found for this genome</em></p>'

        # Build KCB (KnownClusterBlast) section if tabulation data available
        kcb_section = ''
        if tab_df is not None and 'KCB_hit' in tab_df.columns:
            # Filter tabulation data for this genome
            genome_tab = tab_df[tab_df['file'] == genome_name]
            total_regions = len(genome_tab)
            kcb_hits = genome_tab[genome_tab['KCB_hit'] != '']
            regions_with_hits = len(kcb_hits)
            regions_without_hits = total_regions - regions_with_hits
            unique_clusters = kcb_hits['KCB_hit'].nunique() if not kcb_hits.empty else 0

            kcb_section = f'''
    <h2>Known Cluster Matches (KnownClusterBlast)</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>Total BGC Regions</td><td>{total_regions}</td></tr>
        <tr><td>Regions with Known Cluster Matches</td><td>{regions_with_hits}</td></tr>
        <tr><td>Regions without Known Cluster Matches (Potentially Novel)</td><td>{regions_without_hits}</td></tr>
        <tr><td>Unique Known Clusters Matched</td><td>{unique_clusters}</td></tr>
    </table>
    <p style="color: #666; margin-top: 15px;"><em>View detailed cluster matches in the antiSMASH results for this genome.</em></p>
    '''

        html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>Genome Metadata - {genome_name}</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #f8f9fa; color: #333; }}
        h1 {{ color: #333; }}
        h2 {{ color: #333; margin-top: 30px; border-bottom: 2px solid #5b8ac5; padding-bottom: 10px; }}
        table {{ border-collapse: collapse; width: 100%; background: white; margin-top: 20px; border: 1px solid #ddd; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background: #2c5aa0; color: white; font-weight: bold; }}
        tr:nth-child(even) {{ background-color: #f8f9fa; }}
        tr:hover {{ background-color: #e9ecef; }}
        .back-link {{ margin-top: 20px; display: inline-block; padding: 10px 20px; background: #2c5aa0; color: white; font-weight: bold; text-decoration: none; border-radius: 4px; }}
        .back-link:hover {{ background: #5b8ac5; }}
        .antismash-link {{ display: inline-block; padding: 10px 20px; background: #28a745; color: white; font-weight: bold; text-decoration: none; border-radius: 4px; margin-right: 10px; }}
        .antismash-link:hover {{ background: #34ce57; color: white; text-decoration: none; }}
        .button-container {{ margin-top: 30px; }}
        h3 {{ color: #333; margin-top: 25px; }}
        table a {{ color: #2c5aa0; text-decoration: none; }}
        table a:hover {{ text-decoration: underline; }}
    </style>
</head>
<body>
    <h1>Genome Metadata: {genome_name}</h1>

    <h2>BGC Statistics</h2>
    <table>
        <tr><th>Field</th><th>Value</th></tr>
        <tr><td>Total BGC Count</td><td>{row.get('total_count', 'N/A')}</td></tr>
    </table>

    <h2>Assembly Metadata</h2>
    {assembly_table}

    {kcb_section}

    <div class="button-container">
        <a href="../../../../antismash_results/{taxon_clean}/{genome_name}/index.html" class="antismash-link" target="_blank">View antiSMASH Results</a>
        <a href="../bgc_report.html" class="back-link">← Back to Main Report</a>
    </div>
</body>
</html>
'''

        output_file = genome_dir / f'{genome_name}.html'
        with open(output_file, 'w') as f:
            f.write(html_content)

    return len(counts_df)


def generate_genome_table_html(counts_file, assembly_info, name_map, taxonomy_map_data=None):
    '''Generate HTML for a searchable genome table'''

    # Read counts to get genome list with BGC data
    counts_df = pd.read_csv(counts_file, sep='\t', comment='#')

    # Read assembly info
    assembly_df = pd.read_csv(assembly_info, sep='\t')

    # Load name map and create reverse mapping (genome_name -> assembly_id)
    with open(name_map, 'r') as f:
        name_map_data = json.load(f)
    reverse_map = {v: k for k, v in name_map_data.items()}

    # Get BGC columns (exclude record and total_count)
    numeric_cols = counts_df.select_dtypes(include='number').columns
    bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

    # Build table rows
    table_rows = []
    for _, row in counts_df.iterrows():
        genome_name = row['record'].replace('.gbff', '')
        assembly_id = reverse_map.get(genome_name, 'N/A')
        total_bgcs = row.get('total_count', 0)

        # Get top BGC types for this genome
        bgc_counts = {col: row.get(col, 0) for col in bgc_cols if row.get(col, 0) > 0}
        top_bgcs = ', '.join([f"{k}({v})" for k, v in sorted(bgc_counts.items(), key=lambda x: -x[1])[:3]])
        if not top_bgcs:
            top_bgcs = '-'

        # Get organism name from assembly info
        metadata_rows = assembly_df[assembly_df.iloc[:, 0] == assembly_id] if assembly_id != 'N/A' else pd.DataFrame()
        organism_name = metadata_rows.iloc[0, 1] if not metadata_rows.empty and len(metadata_rows.columns) > 1 else 'N/A'
        if pd.isna(organism_name):
            organism_name = 'N/A'

        # Get taxonomy info
        taxonomy = 'N/A'
        if taxonomy_map_data and assembly_id in taxonomy_map_data:
            lineage = taxonomy_map_data[assembly_id].get('lineage', {})
            # Get genus and species if available
            genus = lineage.get('genus', {}).get('name', '')
            species = lineage.get('species', {}).get('name', '')
            if species:
                taxonomy = species
            elif genus:
                taxonomy = genus

        table_rows.append({
            'genome_name': genome_name,
            'assembly_id': assembly_id,
            'organism': str(organism_name)[:50],  # Truncate long names
            'taxonomy': taxonomy,
            'total_bgcs': int(total_bgcs),
            'top_bgcs': top_bgcs
        })

    # Sort by total BGCs descending
    table_rows.sort(key=lambda x: -x['total_bgcs'])

    # Generate HTML rows
    html_rows = []
    for r in table_rows:
        html_rows.append(f'''
            <tr>
                <td><a href="genomes/{r['genome_name']}.html">{r['genome_name']}</a></td>
                <td>{r['assembly_id']}</td>
                <td title="{r['organism']}">{r['organism']}</td>
                <td>{r['taxonomy']}</td>
                <td>{r['total_bgcs']}</td>
                <td>{r['top_bgcs']}</td>
            </tr>''')

    return '\n'.join(html_rows)


def calculate_summary_statistics(counts_file, tabulation_file=None):
    '''Calculate summary statistics for the dataset including tabulation stats'''
    df = pd.read_csv(counts_file, sep='\t', comment='#')

    # Select BGC columns
    numeric_cols = df.select_dtypes(include='number').columns
    bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

    total_genomes = len(df)
    total_bgcs = df['total_count'].sum()
    avg_bgcs = df['total_count'].mean()
    std_bgcs = df['total_count'].std()
    min_bgcs = int(df['total_count'].min())
    max_bgcs = int(df['total_count'].max())
    median_bgcs = df['total_count'].median()
    genomes_with_no_bgcs = len(df[df['total_count'] == 0])
    genomes_with_bgcs = total_genomes - genomes_with_no_bgcs

    # Find most common BGC type(s)
    bgc_totals = df[bgc_cols].sum().sort_values(ascending=False)
    most_common_bgc = bgc_totals.index[0] if len(bgc_totals) > 0 else 'N/A'
    most_common_count = bgc_totals.iloc[0] if len(bgc_totals) > 0 else 0

    # Calculate tabulation stats if file provided
    kcb_stats = {
        'total_regions': 0,
        'regions_with_hits': 0,
        'hit_percentage': 0,
        'unique_known_clusters': 0,
        'similarity_breakdown': {},
        'contig_edge_count': 0,
        'top_known_clusters': []
    }

    if tabulation_file and Path(tabulation_file).exists():
        try:
            tab_df = pd.read_csv(tabulation_file, sep='\t')
            kcb_stats['total_regions'] = len(tab_df)

            # Fill NaN values with empty strings for consistent filtering
            for col in ['KCB_hit', 'KCB_acc', 'KCB_sim']:
                if col in tab_df.columns:
                    tab_df[col] = tab_df[col].fillna('')

            # Count regions with KnownClusterBlast hits
            if 'KCB_sim' in tab_df.columns:
                regions_with_hits = len(tab_df[tab_df['KCB_sim'] != ''])
                kcb_stats['regions_with_hits'] = int(regions_with_hits)
                kcb_stats['hit_percentage'] = round(regions_with_hits / len(tab_df) * 100, 1) if len(tab_df) > 0 else 0

                # Similarity breakdown
                sim_counts = tab_df[tab_df['KCB_sim'] != '']['KCB_sim'].value_counts().to_dict()
                kcb_stats['similarity_breakdown'] = sim_counts

            # Count unique known clusters
            if 'KCB_hit' in tab_df.columns:
                unique_clusters = tab_df[tab_df['KCB_hit'] != '']['KCB_hit'].nunique()
                kcb_stats['unique_known_clusters'] = int(unique_clusters)

            # Build known cluster to region mapping
            if 'KCB_hit' in tab_df.columns and 'KCB_acc' in tab_df.columns:
                kcb_mapping = []
                hits_df = tab_df[tab_df['KCB_hit'] != ''].copy()
                if len(hits_df) > 0:
                    # Group by known cluster
                    for kcb_name in hits_df['KCB_hit'].unique():
                        cluster_rows = hits_df[hits_df['KCB_hit'] == kcb_name]
                        regions = []
                        for _, row in cluster_rows.iterrows():
                            genome = row.get('file', 'unknown')
                            region_num = row.get('region', '?')
                            region_name = row.get('region_name', region_num)
                            record_index = row.get('record_index', 1)
                            product = row.get('product', '')
                            sim = row.get('KCB_sim', '')
                            regions.append({
                                'genome': genome,
                                'region': region_num,
                                'region_name': region_name,
                                'record_index': record_index,
                                'product': product,
                                'similarity': sim
                            })
                        kcb_acc = cluster_rows.iloc[0]['KCB_acc'] if 'KCB_acc' in cluster_rows.columns else ''
                        kcb_mapping.append({
                            'known_cluster': kcb_name,
                            'mibig_acc': kcb_acc,
                            'regions': regions,
                            'count': len(regions)
                        })
                    # Sort by count descending
                    kcb_mapping.sort(key=lambda x: x['count'], reverse=True)
                kcb_stats['cluster_mapping'] = kcb_mapping

            # Build novel BGCs list (regions without KCB hits)
            if 'KCB_hit' in tab_df.columns:
                novel_df = tab_df[tab_df['KCB_hit'] == ''].copy()
                novel_bgcs = []
                if len(novel_df) > 0:
                    for _, row in novel_df.iterrows():
                        novel_bgcs.append({
                            'genome': row.get('file', 'unknown'),
                            'region': row.get('region', '?'),
                            'region_name': row.get('region_name', row.get('region', '?')),
                            'record_index': row.get('record_index', 1),
                            'product': row.get('product', ''),
                            'record_id': row.get('record_id', ''),
                            'contig_edge': row.get('contig_edge', '')
                        })
                kcb_stats['novel_bgcs'] = novel_bgcs
                kcb_stats['novel_bgc_count'] = len(novel_bgcs)

            # Count BGCs on contig edges (potentially incomplete)
            if 'contig_edge' in tab_df.columns:
                contig_edge_count = len(tab_df[tab_df['contig_edge'].astype(str).str.lower() == 'true'])
                kcb_stats['contig_edge_count'] = int(contig_edge_count)

            # Get top 5 known clusters by frequency
            if 'cluster_mapping' in kcb_stats and kcb_stats['cluster_mapping']:
                kcb_stats['top_known_clusters'] = kcb_stats['cluster_mapping'][:5]

        except Exception as e:
            print(f"Warning: Could not parse tabulation file for stats: {e}")

    return {
        'total_genomes': total_genomes,
        'total_bgcs': int(total_bgcs),
        'avg_bgcs': f'{avg_bgcs:.1f}',
        'std_bgcs': f'{std_bgcs:.1f}' if not pd.isna(std_bgcs) else '0.0',
        'min_bgcs': min_bgcs,
        'max_bgcs': max_bgcs,
        'median_bgcs': f'{median_bgcs:.1f}',
        'genomes_with_no_bgcs': genomes_with_no_bgcs,
        'genomes_with_bgcs': genomes_with_bgcs,
        'most_common_bgc': most_common_bgc,
        'most_common_count': int(most_common_count),
        'kcb_stats': kcb_stats
    }

def create_bgc_distribution_table(counts_file, outdir):
    '''Create interactive HTML table for BGC distribution with clickable genome links and color-coding'''
    df = pd.read_csv(counts_file, sep='\t', comment='#')

    # Select BGC columns
    numeric_cols = df.select_dtypes(include='number').columns
    bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

    # Calculate max values for color scaling
    max_total = df['total_count'].max() if df['total_count'].max() > 0 else 1
    max_values = {}
    for col in bgc_cols:
        max_values[col] = df[col].max() if df[col].max() > 0 else 1

    def get_color(value, max_val):
        '''Generate color based on value using sequential green palette'''
        if value == 0 or pd.isna(value):
            return ''
        # Scale from light green to dark green
        intensity = value / max_val
        # Sequential green palette: light (#e8f5e9) to dark (#1b5e20)
        # Light green: rgb(232, 245, 233) Dark green: rgb(27, 94, 32)
        r = int(232 - (205 * intensity))  # 232 to 27
        g = int(245 - (151 * intensity))  # 245 to 94
        b = int(233 - (201 * intensity))  # 233 to 32
        text_color = '#fff' if intensity > 0.6 else '#333'
        return f'background-color: rgb({r}, {g}, {b}); color: {text_color}; font-weight: bold;'

    # Create HTML table
    html_rows = []
    for _, row in df.iterrows():
        genome_name = row['record'].replace('.gbff', '')
        genome_link = f'genomes/{genome_name}.html'

        # Create row with clickable genome name
        row_html = f'<tr><td><a href="{genome_link}">{genome_name}</a></td>'

        # Total count with color
        total_color = get_color(row["total_count"], max_total)
        row_html += f'<td style="{total_color}">{int(row["total_count"])}</td>'

        # Add BGC counts for each type with color coding
        for bgc_type in bgc_cols:
            count = row.get(bgc_type, 0)
            color = get_color(count, max_values[bgc_type])
            cell_value = int(count) if count > 0 else ''
            row_html += f'<td style="{color}">{cell_value}</td>'

        row_html += '</tr>'
        html_rows.append(row_html)

    # Create header
    header = '<tr><th>Genome (click for details)</th><th>Total BGCs</th>'
    for bgc_type in bgc_cols:
        header += f'<th>{bgc_type}</th>'
    header += '</tr>'

    return header, '\n'.join(html_rows)


def plot_bgc_donut_chart(counts_file, outdir):
    '''Create horizontal bar chart showing overall BGC type distribution with consistent colors'''
    df = pd.read_csv(counts_file, sep='\t', comment='#')

    # Select BGC columns
    numeric_cols = df.select_dtypes(include='number').columns
    bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

    # Sum total counts for each BGC type
    bgc_totals = df[bgc_cols].sum().sort_values(ascending=False)

    # Filter out zero values
    bgc_totals = bgc_totals[bgc_totals > 0]

    if len(bgc_totals) == 0:
        print("No BGC data for distribution chart")
        return False

    # Calculate percentages and filter out entries below 0.5% (would round to 0%)
    total_bgcs = bgc_totals.sum()
    bgc_percentages = (bgc_totals / total_bgcs * 100)

    # Filter to only include entries with >= 0.5% (visible percentages)
    significant_mask = bgc_percentages >= 0.5
    bgc_totals_filtered = bgc_totals[significant_mask]
    bgc_percentages_filtered = bgc_percentages[significant_mask]

    # Track filtered entries for note
    n_filtered = len(bgc_totals) - len(bgc_totals_filtered)
    filtered_sum = bgc_totals[~significant_mask].sum() if n_filtered > 0 else 0

    if len(bgc_totals_filtered) == 0:
        print("No BGC data above 0.5% threshold for distribution chart")
        return False

    # Get colors for each BGC type
    colors = [get_bgc_color(bgc_type) for bgc_type in bgc_totals_filtered.index]

    # Create horizontal bar chart - wider figure for full-width display
    n_types = len(bgc_totals_filtered)
    fig_height = max(6, n_types * 0.4 + 2)  # Dynamic height based on number of types
    fig, ax = plt.subplots(figsize=(14, fig_height), facecolor='white')
    ax.set_facecolor('white')

    # Create horizontal bars
    y_pos = range(len(bgc_totals_filtered))
    bars = ax.barh(y_pos, bgc_percentages_filtered.values, color=colors, edgecolor='white', linewidth=1)

    # Add value labels on bars
    for i, (bar, count, pct) in enumerate(zip(bars, bgc_totals_filtered.values, bgc_percentages_filtered.values)):
        # Label inside bar if wide enough, otherwise outside
        if pct > 5:
            ax.text(bar.get_width() - 0.5, bar.get_y() + bar.get_height()/2,
                   f'{int(count)} ({pct:.1f}%)', ha='right', va='center',
                   fontsize=11, fontweight='bold', color='white')
        else:
            ax.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                   f'{int(count)} ({pct:.1f}%)', ha='left', va='center',
                   fontsize=11, fontweight='bold', color='#333')

    # Styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels(bgc_totals_filtered.index, fontsize=12)
    ax.set_xlabel('Percentage of Total BGCs', fontsize=12, fontweight='bold')
    ax.set_xlim(0, max(bgc_percentages_filtered.values) * 1.15)  # Add space for labels
    ax.invert_yaxis()  # Largest at top

    # Add title with total count
    title = f'BGC Type Distribution (n={int(total_bgcs)} total BGCs)'
    if n_filtered > 0:
        title += f'\n{n_filtered} types with <0.5% ({int(filtered_sum)} BGCs) not shown'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=10)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{outdir}/bgc_donut_chart.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return True

def plot_circular_taxonomy_tree(taxonomy_tree_data, counts_file, outdir):
    '''Create circular taxonomy tree with BGC annotation rings'''
    if not taxonomy_tree_data:
        print("No taxonomy tree data for circular visualization")
        return False

    tree = taxonomy_tree_data.get('tree', {})
    if not tree:
        print("Empty taxonomy tree")
        return False

    # Read counts data to get BGC types
    df = pd.read_csv(counts_file, sep='\t', comment='#')
    numeric_cols = df.select_dtypes(include='number').columns
    bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

    # Collect all leaf nodes (genomes) with their taxonomic path
    leaves = []
    taxonomic_groups = defaultdict(list)  # group_key -> list of leaf indices

    def collect_leaves(node, path=None, level=0):
        '''Recursively collect leaf nodes with their taxonomic paths'''
        if path is None:
            path = []

        name = node.get('name', 'Unknown')
        rank = node.get('rank', '')
        children = node.get('children', {})
        genomes = node.get('genomes', [])

        current_path = path + [(name, rank, level)]

        # If this is a species node with genomes, add them as leaves
        if genomes and rank == 'species':
            for genome in genomes:
                genome_name = genome.get('name', 'Unknown')
                total_bgcs = genome.get('total_bgcs', 0)
                bgc_types = genome.get('bgc_types', {})
                leaf_idx = len(leaves)
                leaves.append({
                    'name': genome_name,
                    'path': current_path,
                    'total_bgcs': total_bgcs,
                    'bgc_types': bgc_types
                })
                # Track taxonomic groupings for arc drawing
                for i, (tax_name, tax_rank, tax_level) in enumerate(current_path):
                    group_key = (tax_name, tax_rank, tax_level)
                    taxonomic_groups[group_key].append(leaf_idx)

        # Recurse into children
        for child_name, child_node in children.items():
            collect_leaves(child_node, current_path, level + 1)

    collect_leaves(tree)

    if len(leaves) == 0:
        print("No leaf nodes found in taxonomy tree")
        return False

    print(f"Found {len(leaves)} genomes for circular tree")

    # Calculate angular positions for leaves
    n_leaves = len(leaves)
    angles = np.linspace(0, 2 * np.pi, n_leaves, endpoint=False)

    # Assign angles to leaves
    for i, leaf in enumerate(leaves):
        leaf['angle'] = angles[i]

    # Create figure with polar projection
    fig = plt.figure(figsize=(16, 16), facecolor='white')

    # Main circular tree axis
    ax = fig.add_subplot(111, projection='polar')
    ax.set_facecolor('white')

    # Tree parameters
    max_level = max(len(leaf['path']) for leaf in leaves)
    inner_radius = 0.3
    tree_radius = 0.6
    level_step = (tree_radius - inner_radius) / (max_level + 1)

    # Draw taxonomic arcs (from innermost to outermost)
    # Group leaves by taxonomic level and draw arcs
    arc_colors = {
        'domain': '#e8e8e8',
        'kingdom': '#d0d0d0',
        'phylum': '#b8b8b8',
        'class': '#a0a0a0',
        'order': '#888888',
        'family': '#707070',
        'genus': '#585858',
        'species': '#404040'
    }

    for (tax_name, tax_rank, tax_level), leaf_indices in taxonomic_groups.items():
        if len(leaf_indices) < 1:
            continue

        # Get angles for this group
        group_angles = [leaves[i]['angle'] for i in leaf_indices]
        min_angle = min(group_angles)
        max_angle = max(group_angles)

        # Handle wrap-around case
        angle_span = max_angle - min_angle
        if angle_span > np.pi:
            # Wrap around - need to handle differently
            angles_sorted = sorted(group_angles)
            gaps = [angles_sorted[i+1] - angles_sorted[i] for i in range(len(angles_sorted)-1)]
            gaps.append(2*np.pi - angles_sorted[-1] + angles_sorted[0])
            max_gap_idx = gaps.index(max(gaps))
            if max_gap_idx == len(gaps) - 1:
                min_angle = angles_sorted[0]
                max_angle = angles_sorted[-1]
            else:
                min_angle = angles_sorted[max_gap_idx + 1]
                max_angle = angles_sorted[max_gap_idx]
                if max_angle < min_angle:
                    max_angle += 2 * np.pi

        # Calculate radius for this level
        radius = inner_radius + (tax_level + 0.5) * level_step

        # Draw arc
        arc_color = arc_colors.get(tax_rank, '#cccccc')

        # Add small padding to angles
        pad = 0.01
        theta1 = min_angle - pad
        theta2 = max_angle + pad

        # Draw arc segment
        arc_angles = np.linspace(theta1, theta2, 50)
        ax.plot(arc_angles, [radius] * len(arc_angles), color=arc_color, linewidth=3, alpha=0.7)

    # Draw radial lines connecting leaves to center
    for leaf in leaves:
        angle = leaf['angle']
        ax.plot([angle, angle], [inner_radius, tree_radius], color='#e0e0e0', linewidth=0.5, alpha=0.5)

    # Draw BGC count ring (outer ring showing total BGCs per genome)
    bgc_ring_inner = tree_radius + 0.05
    bgc_ring_outer = tree_radius + 0.15

    # Normalize BGC counts for ring height
    max_bgcs = max(leaf['total_bgcs'] for leaf in leaves) if leaves else 1
    max_bgcs = max(max_bgcs, 1)

    # Draw BGC bars for each genome
    bar_width = 2 * np.pi / n_leaves * 0.8

    for leaf in leaves:
        angle = leaf['angle']
        total = leaf['total_bgcs']

        if total > 0:
            # Calculate bar height based on BGC count
            bar_height = (total / max_bgcs) * (bgc_ring_outer - bgc_ring_inner)

            # Draw stacked bar by BGC type
            current_bottom = bgc_ring_inner
            for bgc_type in bgc_cols:
                count = leaf['bgc_types'].get(bgc_type, 0)
                if count > 0:
                    segment_height = (count / total) * bar_height
                    color = get_bgc_color(bgc_type)

                    # Draw wedge for this BGC type
                    ax.bar(angle, segment_height, width=bar_width, bottom=current_bottom,
                           color=color, edgecolor='white', linewidth=0.2, alpha=0.9)
                    current_bottom += segment_height

    # Draw genome labels (only if not too many)
    if n_leaves <= 100:
        label_radius = bgc_ring_outer + 0.08
        for leaf in leaves:
            angle = leaf['angle']
            name = leaf['name']

            # Truncate long names
            if len(name) > 25:
                name = name[:22] + '...'

            # Rotate text to be readable
            rotation = np.degrees(angle) - 90
            if angle > np.pi/2 and angle < 3*np.pi/2:
                rotation += 180
                ha = 'right'
            else:
                ha = 'left'

            ax.text(angle, label_radius, name, rotation=rotation, ha=ha, va='center',
                    fontsize=6, color='#333', rotation_mode='anchor')

    # Configure polar plot
    ax.set_ylim(0, bgc_ring_outer + 0.3)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.axis('off')

    # Add title
    fig.suptitle('Circular Taxonomy Tree with BGC Distribution',
                 fontsize=18, fontweight='bold', color='#333', y=0.98)

    # Add legend for BGC types (only show types present in data)
    present_bgc_types = set()
    for leaf in leaves:
        present_bgc_types.update(leaf['bgc_types'].keys())

    legend_handles = []
    for bgc_type in sorted(present_bgc_types):
        color = get_bgc_color(bgc_type)
        patch = mpatches.Patch(color=color, label=bgc_type)
        legend_handles.append(patch)

    if legend_handles:
        # Position legend outside the plot
        legend = fig.legend(handles=legend_handles, title='BGC Types',
                           loc='center left', bbox_to_anchor=(0.85, 0.5),
                           fontsize=8, title_fontsize=10)
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_edgecolor('#ddd')

    # Add note about what the rings represent
    fig.text(0.02, 0.02,
             'Inner rings: Taxonomic hierarchy (darker = lower rank)\n'
             'Outer ring: BGC counts per genome (stacked by type, height = total BGCs)',
             fontsize=9, color='#666', va='bottom', ha='left',
             transform=fig.transFigure)

    plt.tight_layout(rect=[0, 0.03, 0.85, 0.97])
    plt.savefig(f'{outdir}/circular_taxonomy_tree.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


# =============================================================================
# PHYLOGENETIC TREE VISUALIZATION (From Newick format - GTDB-Tk output)
# =============================================================================

def parse_newick(newick_str):
    """
    Parse a Newick format string into a tree structure.
    Returns a dictionary representing the tree with nodes containing:
    - name: node name (leaf names or internal node labels)
    - branch_length: distance from parent
    - children: list of child nodes
    """
    import re

    # Clean up the string
    newick_str = newick_str.strip()
    if newick_str.endswith(';'):
        newick_str = newick_str[:-1]

    def parse_subtree(s, pos=0):
        """Recursively parse a subtree from position pos"""
        node = {'name': '', 'branch_length': 0.0, 'children': []}

        if pos >= len(s):
            return node, pos

        # Check if this is an internal node (starts with '(')
        if s[pos] == '(':
            pos += 1  # skip '('

            # Parse children
            while True:
                child, pos = parse_subtree(s, pos)
                node['children'].append(child)

                if pos >= len(s):
                    break

                if s[pos] == ',':
                    pos += 1  # skip ',' and continue to next child
                elif s[pos] == ')':
                    pos += 1  # skip ')' and break
                    break

        # Parse node name (can be after ')' for internal nodes or at current position for leaves)
        name_match = re.match(r"([^:,();\[\]]*)", s[pos:])
        if name_match:
            node['name'] = name_match.group(1).strip().strip("'\"")
            pos += len(name_match.group(0))

        # Parse branch length if present
        if pos < len(s) and s[pos] == ':':
            pos += 1
            length_match = re.match(r"([0-9.eE+-]+)", s[pos:])
            if length_match:
                try:
                    node['branch_length'] = float(length_match.group(1))
                except ValueError:
                    node['branch_length'] = 0.0
                pos += len(length_match.group(0))

        # Skip any trailing metadata in brackets (e.g., bootstrap values)
        if pos < len(s) and s[pos] == '[':
            bracket_count = 1
            pos += 1
            while pos < len(s) and bracket_count > 0:
                if s[pos] == '[':
                    bracket_count += 1
                elif s[pos] == ']':
                    bracket_count -= 1
                pos += 1

        return node, pos

    tree, _ = parse_subtree(newick_str)
    return tree


def collect_phylo_leaves(node, depth=0.0):
    """
    Recursively collect leaf nodes with their cumulative distances from root.
    Returns list of (leaf_name, cumulative_distance) tuples.
    """
    leaves = []
    current_depth = depth + node.get('branch_length', 0.0)

    if not node.get('children'):
        # This is a leaf node
        return [(node.get('name', ''), current_depth, node)]

    for child in node.get('children', []):
        leaves.extend(collect_phylo_leaves(child, current_depth))

    return leaves


def prune_tree_to_leaves(node, target_leaves):
    """
    Prune a tree to only include paths from root to specified leaf nodes.

    Args:
        node: Tree node dictionary with 'name', 'branch_length', 'children'
        target_leaves: Set of leaf names to keep

    Returns:
        Pruned tree node, or None if this subtree contains no target leaves
    """
    # If this is a leaf node
    if not node.get('children'):
        leaf_name = node.get('name', '')
        # Check if this leaf should be kept (exact match or partial match)
        for target in target_leaves:
            if leaf_name == target or target in leaf_name or leaf_name in target:
                return {
                    'name': node.get('name', ''),
                    'branch_length': node.get('branch_length', 0.0),
                    'children': []
                }
        return None

    # Internal node - recursively prune children
    pruned_children = []
    for child in node.get('children', []):
        pruned_child = prune_tree_to_leaves(child, target_leaves)
        if pruned_child is not None:
            pruned_children.append(pruned_child)

    # If no children remain after pruning, this subtree is not needed
    if not pruned_children:
        return None

    # If only one child remains, we can optionally collapse this node
    # But for phylogenetic trees, we want to preserve branch lengths
    # so we keep the node structure

    return {
        'name': node.get('name', ''),
        'branch_length': node.get('branch_length', 0.0),
        'children': pruned_children
    }


def tree_to_newick(node):
    """
    Convert a tree dictionary back to Newick format string.

    Args:
        node: Tree node dictionary with 'name', 'branch_length', 'children'

    Returns:
        Newick format string
    """
    if not node:
        return ""

    children = node.get('children', [])
    name = node.get('name', '')
    branch_length = node.get('branch_length', 0.0)

    if not children:
        # Leaf node
        if branch_length > 0:
            return f"{name}:{branch_length}"
        return name

    # Internal node
    child_strings = [tree_to_newick(child) for child in children]
    children_str = ','.join(child_strings)

    if branch_length > 0:
        return f"({children_str}){name}:{branch_length}"
    elif name:
        return f"({children_str}){name}"
    else:
        return f"({children_str})"


def calculate_phylo_positions(node, leaf_angles, depth=0.0, angle_range=None):
    """
    Recursively calculate angular positions for all nodes.
    Returns dict mapping node names to (angle, radius) tuples.
    """
    positions = {}
    current_depth = depth + node.get('branch_length', 0.0)

    if not node.get('children'):
        # Leaf node - get pre-assigned angle
        name = node.get('name', '')
        if name in leaf_angles:
            positions[name] = (leaf_angles[name], current_depth)
        return positions

    # Internal node - calculate angle as average of children's angles
    child_angles = []
    for child in node.get('children', []):
        child_positions = calculate_phylo_positions(child, leaf_angles, current_depth)
        positions.update(child_positions)

        # Get the child's angle (either directly or from its descendants)
        child_name = child.get('name', '')
        if child_name in positions:
            child_angles.append(positions[child_name][0])
        elif child.get('children'):
            # Average of grandchildren angles
            grandchild_angles = [positions[gc.get('name', '')][0]
                                for gc in child.get('children', [])
                                if gc.get('name', '') in positions]
            if grandchild_angles:
                child_angles.append(np.mean(grandchild_angles))

    # Position internal node at average of children angles
    if child_angles:
        node_angle = np.mean(child_angles)
        node_name = node.get('name', '') or f"_internal_{id(node)}"
        positions[node_name] = (node_angle, current_depth)

    return positions


def plot_circular_phylogenetic_tree(newick_file, counts_file, gtdbtk_summary, outdir):
    """
    Create circular phylogenetic tree from GTDB-Tk Newick output with BGC annotation rings.

    Args:
        newick_file: Path to Newick tree file from GTDB-Tk
        counts_file: Path to region_counts.tsv for BGC data
        gtdbtk_summary: Path to GTDB-Tk summary TSV for genome name mapping
        outdir: Output directory for the plot
    """
    from pathlib import Path

    newick_path = Path(newick_file)
    if not newick_path.exists():
        print(f"Newick file not found: {newick_file}")
        return False

    # Read GTDB-Tk summary to get user genome names (for pruning if needed)
    user_genomes = set()
    genome_name_map = {}  # Maps GTDB-Tk IDs to our genome names
    if gtdbtk_summary and Path(gtdbtk_summary).exists():
        try:
            gtdbtk_df = pd.read_csv(gtdbtk_summary, sep='\t')
            if 'user_genome' in gtdbtk_df.columns:
                for _, row in gtdbtk_df.iterrows():
                    user_genome = row['user_genome']
                    user_genomes.add(user_genome)
                    genome_name_map[user_genome] = user_genome
            print(f"Found {len(user_genomes)} user genomes in GTDB-Tk summary")
        except Exception as e:
            print(f"Warning: Could not read GTDB-Tk summary: {e}")

    # Read and parse Newick tree
    print(f"Parsing phylogenetic tree from {newick_file}...")
    with open(newick_file, 'r') as f:
        newick_str = f.read()

    tree = parse_newick(newick_str)

    if not tree:
        print("Failed to parse Newick tree")
        return False

    # Count leaves in tree
    leaves_data = collect_phylo_leaves(tree)
    n_leaves = len(leaves_data)
    print(f"Tree has {n_leaves} leaves")

    # Only prune if tree has significantly more leaves than user genomes
    # (de_novo_wf with --skip_gtdb_refs produces trees with only user genomes)
    if user_genomes and n_leaves > len(user_genomes) * 1.5:
        print(f"Tree has more leaves than user genomes - pruning to {len(user_genomes)} user genomes...")
        tree = prune_tree_to_leaves(tree, user_genomes)

        if not tree:
            print("Failed to prune tree - no user genomes found in tree")
            return False

        # Recollect leaves after pruning
        leaves_data = collect_phylo_leaves(tree)
        n_leaves = len(leaves_data)
        print(f"Pruned tree has {n_leaves} leaves")

    if n_leaves == 0:
        print("No leaf nodes found in phylogenetic tree")
        return False

    # Read BGC counts data
    bgc_data = {}
    bgc_cols = []
    if counts_file and Path(counts_file).exists():
        df = pd.read_csv(counts_file, sep='\t', comment='#')
        numeric_cols = df.select_dtypes(include='number').columns
        bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

        for _, row in df.iterrows():
            genome_name = row.get('genome', row.get('file', ''))
            if genome_name:
                # Clean up genome name for matching
                clean_name = Path(genome_name).stem if '/' in str(genome_name) else genome_name
                clean_name = clean_name.replace('.gbff', '').replace('.fna', '')
                bgc_data[clean_name] = {
                    'total_bgcs': row.get('total_count', 0),
                    'bgc_types': {col: row.get(col, 0) for col in bgc_cols if row.get(col, 0) > 0}
                }

    # Assign angular positions to leaves (evenly spaced)
    leaf_angles = {}
    for i, (leaf_name, depth, node) in enumerate(leaves_data):
        angle = 2 * np.pi * i / n_leaves
        # Clean leaf name for matching
        clean_name = leaf_name.replace('.fna', '').replace('.gbff', '')
        leaf_angles[leaf_name] = angle
        if clean_name not in bgc_data:
            # Try to find matching genome
            for genome_key in bgc_data.keys():
                if clean_name in genome_key or genome_key in clean_name:
                    bgc_data[clean_name] = bgc_data[genome_key]
                    break

    # Calculate positions for all nodes
    positions = calculate_phylo_positions(tree, leaf_angles)

    # Find maximum depth for scaling
    max_depth = max(depth for _, depth, _ in leaves_data) if leaves_data else 1.0
    max_depth = max(max_depth, 0.001)  # Avoid division by zero

    # Create figure
    fig = plt.figure(figsize=(18, 18), facecolor='white')
    ax = fig.add_subplot(111, projection='polar')
    ax.set_facecolor('white')

    # Tree layout parameters
    inner_radius = 0.15  # Start of tree
    tree_radius = 0.55   # End of tree (leaves)
    scale_factor = (tree_radius - inner_radius) / max_depth

    def draw_branch(node, parent_angle=None, parent_radius=None, depth=0.0):
        """Recursively draw tree branches"""
        current_depth = depth + node.get('branch_length', 0.0)
        current_radius = inner_radius + current_depth * scale_factor

        node_name = node.get('name', '')
        if not node_name:
            node_name = f"_internal_{id(node)}"

        # Get this node's angle
        if node_name in positions:
            node_angle = positions[node_name][0]
        elif node_name in leaf_angles:
            node_angle = leaf_angles[node_name]
        else:
            # For internal nodes, average children angles
            child_angles = []
            for child in node.get('children', []):
                child_name = child.get('name', '') or f"_internal_{id(child)}"
                if child_name in positions:
                    child_angles.append(positions[child_name][0])
                elif child_name in leaf_angles:
                    child_angles.append(leaf_angles[child_name])
            node_angle = np.mean(child_angles) if child_angles else 0

        # Draw branch from parent to this node
        if parent_angle is not None and parent_radius is not None:
            # Draw arc at parent's radius from parent angle to this node's angle
            if abs(node_angle - parent_angle) > 0.001:
                # Determine arc direction (shortest path)
                angle_diff = node_angle - parent_angle
                if angle_diff > np.pi:
                    angle_diff -= 2 * np.pi
                elif angle_diff < -np.pi:
                    angle_diff += 2 * np.pi

                arc_angles = np.linspace(parent_angle, parent_angle + angle_diff, 20)
                ax.plot(arc_angles, [parent_radius] * len(arc_angles),
                       color='#505050', linewidth=0.8, alpha=0.8)

            # Draw radial line from parent's arc to current node
            ax.plot([node_angle, node_angle], [parent_radius, current_radius],
                   color='#505050', linewidth=0.8, alpha=0.8)

        # Recurse for children
        for child in node.get('children', []):
            draw_branch(child, node_angle, current_radius, current_depth)

    # Draw the tree
    print("Drawing phylogenetic tree branches...")
    draw_branch(tree)

    # Draw BGC annotation ring
    bgc_ring_inner = tree_radius + 0.05
    bgc_ring_outer = tree_radius + 0.18

    # Calculate max BGCs for scaling
    max_bgcs = max([bgc_data.get(Path(name).stem.replace('.fna', '').replace('.gbff', ''), {}).get('total_bgcs', 0)
                   for name, _, _ in leaves_data]) if leaves_data else 1
    max_bgcs = max(max_bgcs, 1)

    bar_width = 2 * np.pi / n_leaves * 0.85

    print("Drawing BGC annotation ring...")
    present_bgc_types = set()

    for leaf_name, depth, node in leaves_data:
        angle = leaf_angles[leaf_name]
        clean_name = leaf_name.replace('.fna', '').replace('.gbff', '')

        # Try to find BGC data for this leaf
        genome_bgc = bgc_data.get(clean_name, {})
        if not genome_bgc:
            # Try partial matching
            for key, value in bgc_data.items():
                if clean_name in key or key in clean_name:
                    genome_bgc = value
                    break

        total = genome_bgc.get('total_bgcs', 0)

        if total > 0:
            bar_height = (total / max_bgcs) * (bgc_ring_outer - bgc_ring_inner)
            current_bottom = bgc_ring_inner

            for bgc_type, count in genome_bgc.get('bgc_types', {}).items():
                if count > 0:
                    segment_height = (count / total) * bar_height
                    color = get_bgc_color(bgc_type)
                    present_bgc_types.add(bgc_type)

                    ax.bar(angle, segment_height, width=bar_width, bottom=current_bottom,
                          color=color, edgecolor='white', linewidth=0.2, alpha=0.9)
                    current_bottom += segment_height

    # Draw leaf labels (only if not too many)
    if n_leaves <= 80:
        label_radius = bgc_ring_outer + 0.08
        for leaf_name, depth, node in leaves_data:
            angle = leaf_angles[leaf_name]
            display_name = leaf_name.replace('.fna', '').replace('.gbff', '')

            if len(display_name) > 30:
                display_name = display_name[:27] + '...'

            rotation = np.degrees(angle) - 90
            if angle > np.pi/2 and angle < 3*np.pi/2:
                rotation += 180
                ha = 'right'
            else:
                ha = 'left'

            ax.text(angle, label_radius, display_name, rotation=rotation, ha=ha, va='center',
                   fontsize=5, color='#333', rotation_mode='anchor')

    # Configure plot
    ax.set_ylim(0, bgc_ring_outer + 0.35)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.axis('off')

    # Title
    fig.suptitle('Circular Phylogenetic Tree with BGC Distribution\n(GTDB-Tk placement)',
                fontsize=18, fontweight='bold', color='#333', y=0.98)

    # Legend for BGC types
    legend_handles = []
    for bgc_type in sorted(present_bgc_types):
        color = get_bgc_color(bgc_type)
        patch = mpatches.Patch(color=color, label=bgc_type)
        legend_handles.append(patch)

    if legend_handles:
        legend = fig.legend(handles=legend_handles, title='BGC Types',
                           loc='center left', bbox_to_anchor=(0.85, 0.5),
                           fontsize=8, title_fontsize=10)
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_edgecolor('#ddd')

    # Add explanatory note
    fig.text(0.02, 0.02,
            'Branch lengths represent evolutionary distance (GTDB-Tk)\n'
            'Outer ring: BGC counts per genome (stacked by type, height = total BGCs)',
            fontsize=9, color='#666', va='bottom', ha='left',
            transform=fig.transFigure)

    plt.tight_layout(rect=[0, 0.03, 0.85, 0.97])
    plt.savefig(f'{outdir}/circular_phylogenetic_tree.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved circular phylogenetic tree to {outdir}/circular_phylogenetic_tree.png")
    return True


def prepare_phylo_tree_for_js(newick_file, gtdbtk_summary, counts_file, outdir, outgroup=None):
    """
    Prepare phylogenetic tree data for JavaScript visualization.
    Prunes the GTDB-Tk tree to only include user genomes and exports data for JS rendering.
    Uses Bio.Phylo for efficient parsing of large trees.

    Args:
        newick_file: Path to Newick tree file from GTDB-Tk
        gtdbtk_summary: Path to GTDB-Tk summary TSV
        counts_file: Path to region_counts.tsv for BGC data
        outdir: Output directory
        outgroup: Optional outgroup taxon pattern (e.g., "g__Escherichia") to keep in pruned tree

    Returns:
        dict with 'newick' (pruned tree string) and 'metadata' (genome info), or None on failure
    """
    from pathlib import Path
    from io import StringIO

    newick_path = Path(newick_file)
    if not newick_path.exists():
        print(f"Newick file not found: {newick_file}")
        return None

    # Read user genomes from GTDB-Tk summary
    user_genomes = set()
    genome_metadata = {}
    if gtdbtk_summary and Path(gtdbtk_summary).exists():
        try:
            gtdbtk_df = pd.read_csv(gtdbtk_summary, sep='\t')
            if 'user_genome' in gtdbtk_df.columns:
                for _, row in gtdbtk_df.iterrows():
                    user_genome = row['user_genome']
                    user_genomes.add(user_genome)
                    genome_metadata[user_genome] = {
                        'classification': row.get('classification', ''),
                    }
            print(f"Found {len(user_genomes)} user genomes for tree pruning")
        except Exception as e:
            print(f"Warning: Could not read GTDB-Tk summary: {e}")
            return None

    if not user_genomes:
        print("No user genomes found - cannot create tree")
        return None

    # Read BGC counts
    if counts_file and Path(counts_file).exists():
        try:
            df = pd.read_csv(counts_file, sep='\t', comment='#')
            for _, row in df.iterrows():
                genome_name = row.get('genome', row.get('file', ''))
                if genome_name:
                    clean_name = Path(genome_name).stem if '/' in str(genome_name) else genome_name
                    clean_name = clean_name.replace('.gbff', '').replace('.fna', '')
                    if clean_name in genome_metadata:
                        genome_metadata[clean_name]['total_bgcs'] = int(row.get('total_count', 0))
        except Exception as e:
            print(f"Warning: Could not read counts file: {e}")

    # Try to use Bio.Phylo for efficient parsing
    try:
        from Bio import Phylo
        import sys

        # Increase recursion limit for large trees (GTDB reference trees can be very deep)
        old_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(max(old_limit, 15000))

        print(f"Reading phylogenetic tree from {newick_file} using Bio.Phylo...")

        # Parse tree with Bio.Phylo (much faster than custom parser)
        tree = Phylo.read(newick_file, 'newick')

        # Get all terminal names (leaves)
        terminals = tree.get_terminals()
        print(f"Tree has {len(terminals)} leaves")

        # Find user genomes in tree
        user_terminals = [t for t in terminals if t.name in user_genomes]
        print(f"Found {len(user_terminals)} user genomes in tree")

        if not user_terminals:
            print("No user genomes found in tree - check genome name matching")
            # Try partial matching
            for t in terminals[:10]:
                print(f"  Sample terminal: {t.name}")
            return None

        # Find outgroup terminal if specified
        outgroup_terminal = None
        if outgroup:
            # Search for a terminal matching the outgroup pattern
            # Pattern can be taxonomy rank (e.g., "g__Escherichia") or partial name
            for t in terminals:
                if t.name and outgroup in t.name:
                    outgroup_terminal = t
                    print(f"Found outgroup: {t.name} (matching pattern '{outgroup}')")
                    break
            if not outgroup_terminal:
                print(f"Warning: No outgroup found matching pattern '{outgroup}'")

        # Prune tree to user genomes (and outgroup if found)
        terminals_to_keep = len(user_terminals) + (1 if outgroup_terminal else 0)
        print(f"Pruning tree to {terminals_to_keep} terminals...")

        # For small numbers of genomes, build a simple tree with just their relationships
        if len(user_terminals) == 1 and not outgroup_terminal:
            # Single genome without outgroup - create simple tree
            name = user_terminals[0].name
            bl = user_terminals[0].branch_length or 0.0
            pruned_newick = f"({name}:{bl});"
        else:
            # Use Bio.Phylo's built-in pruning - O(n) instead of O(n²) distance matrix
            # Strategy: remove all terminals except user genomes and outgroup
            terminals_to_keep = set(t.name for t in user_terminals)
            if outgroup_terminal:
                terminals_to_keep.add(outgroup_terminal.name)
            non_user_terminals = [t for t in terminals if t.name not in terminals_to_keep]

            print(f"Removing {len(non_user_terminals)} terminals from tree...")

            # Remove non-user terminals in batches with progress reporting
            removed_count = 0
            total_to_remove = len(non_user_terminals)
            report_interval = max(1, total_to_remove // 20)  # Report ~20 times

            for terminal in non_user_terminals:
                try:
                    tree.prune(terminal)
                    removed_count += 1
                    if removed_count % report_interval == 0:
                        print(f"  Pruning progress: {removed_count}/{total_to_remove} ({100*removed_count//total_to_remove}%)")
                except Exception as e:
                    # Terminal may already be removed if it was part of a collapsed branch
                    pass

            print(f"Removed {removed_count} terminals")

            # Collapse single-child internal nodes to clean up the tree
            def collapse_single_child_clades(clade):
                """Recursively collapse internal nodes with single children."""
                if clade.is_terminal():
                    return clade

                # Process children first
                new_clades = []
                for child in clade.clades:
                    collapsed_child = collapse_single_child_clades(child)
                    if collapsed_child is not None:
                        new_clades.append(collapsed_child)

                clade.clades = new_clades

                # If this node has only one child, merge branch lengths
                if len(clade.clades) == 1:
                    child = clade.clades[0]
                    # Add this node's branch length to child
                    if clade.branch_length and child.branch_length:
                        child.branch_length += clade.branch_length
                    elif clade.branch_length:
                        child.branch_length = clade.branch_length
                    return child

                # If no children left, return None
                if len(clade.clades) == 0:
                    return None

                return clade

            print("Collapsing single-child internal nodes...")
            tree.root = collapse_single_child_clades(tree.root)

            # Write to Newick
            output = StringIO()
            Phylo.write(tree, output, 'newick')
            pruned_newick = output.getvalue().strip()

        print(f"Pruned tree newick length: {len(pruned_newick)} chars")

        # Save pruned Newick to file
        pruned_newick_path = Path(outdir) / 'pruned_phylo_tree.nwk'
        with open(pruned_newick_path, 'w') as f:
            f.write(pruned_newick)
        print(f"Saved pruned Newick to {pruned_newick_path}")

        # Count leaves in pruned tree
        pruned_tree = Phylo.read(StringIO(pruned_newick), 'newick')
        pruned_leaf_count = len(pruned_tree.get_terminals())

        return {
            'newick': pruned_newick,
            'metadata': genome_metadata,
            'leaf_count': pruned_leaf_count
        }

    except ImportError:
        print("Bio.Phylo not available, falling back to custom parser...")
    except Exception as e:
        print(f"Bio.Phylo parsing failed: {e}, falling back to custom parser...")

    # Fallback to custom parser for small trees
    print(f"Reading phylogenetic tree from {newick_file}...")
    with open(newick_file, 'r') as f:
        newick_str = f.read()

    # Quick estimate of tree size
    import re
    leaf_pattern = re.compile(r'[(),]([A-Za-z_][^:(),]*):')
    estimated_leaves = len(leaf_pattern.findall(newick_str))
    print(f"Estimated tree size: ~{estimated_leaves} leaves")

    # Parse Newick tree with custom parser
    print(f"Parsing phylogenetic tree...")
    tree = parse_newick(newick_str)
    if not tree:
        print("Failed to parse Newick tree")
        return None

    original_leaves = collect_phylo_leaves(tree)
    print(f"Original tree has {len(original_leaves)} leaves")

    print(f"Pruning tree to {len(user_genomes)} user genomes...")
    pruned_tree = prune_tree_to_leaves(tree, user_genomes)

    if not pruned_tree:
        print("Failed to prune tree - no user genomes found")
        return None

    pruned_leaves = collect_phylo_leaves(pruned_tree)
    print(f"Pruned tree has {len(pruned_leaves)} leaves")

    pruned_newick = tree_to_newick(pruned_tree) + ";"

    pruned_newick_path = Path(outdir) / 'pruned_phylo_tree.nwk'
    with open(pruned_newick_path, 'w') as f:
        f.write(pruned_newick)
    print(f"Saved pruned Newick to {pruned_newick_path}")

    return {
        'newick': pruned_newick,
        'metadata': genome_metadata,
        'leaf_count': len(pruned_leaves)
    }


def plot_kcb_identification_chart(kcb_stats, outdir):
    '''Create pie chart showing proportion of identified vs unidentified BGCs'''
    total_regions = kcb_stats.get('total_regions', 0)
    regions_with_hits = kcb_stats.get('regions_with_hits', 0)

    if total_regions == 0:
        # No data available - create placeholder (same size as donut chart)
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
        ax.set_facecolor('white')
        ax.text(0.5, 0.5, 'No KnownClusterBlast data available\nRun with --antismash_cb_knownclusters true',
                ha='center', va='center', fontsize=16, color='#666', transform=ax.transAxes)
        ax.axis('off')
        plt.savefig(f'{outdir}/kcb_identification_chart.png', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        return

    regions_without_hits = total_regions - regions_with_hits

    # Same size as donut chart for consistency
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
    ax.set_facecolor('white')

    # Create pie chart with identified vs unidentified
    sizes = [regions_with_hits, regions_without_hits]
    labels = [f'Known Cluster Matches\n({regions_with_hits})', f'Potentially Novel\n({regions_without_hits})']
    colors = ['#2c5aa0', '#4a8f70']  # Deep blue for known, muted green for novel
    explode = (0.02, 0.02)  # Slight separation

    wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                                       startangle=90, explode=explode,
                                       textprops={'fontsize': 14, 'color': '#333'},
                                       wedgeprops={'edgecolor': 'white', 'linewidth': 2})

    # Style the percentage labels
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(16)

    # Add similarity breakdown as legend if available
    sim_breakdown = kcb_stats.get('similarity_breakdown', {})
    if sim_breakdown:
        legend_text = 'Similarity breakdown of matches:\n'
        for sim_level in ['high', 'medium', 'low']:
            if sim_level in sim_breakdown:
                legend_text += f'  • {sim_level.capitalize()}: {sim_breakdown[sim_level]}\n'
        ax.annotate(legend_text.strip(), xy=(0.5, -0.12), xycoords='axes fraction',
                    ha='center', va='top', fontsize=12, color='#666',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8f9fa', edgecolor='#ddd'))

    plt.tight_layout()
    plt.savefig(f'{outdir}/kcb_identification_chart.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def generate_rarefaction_curve(bigscape_db_path, outdir, taxon, n_iterations=50):
    """
    Generate GCF rarefaction curve from BiG-SCAPE database.

    Shows how the number of unique Gene Cluster Families (GCFs) discovered
    increases as more genomes are sampled.

    Returns:
        dict with rarefaction statistics, or None if generation failed
    """
    import os
    import sqlite3
    from collections import defaultdict
    import random

    if not bigscape_db_path or not os.path.exists(bigscape_db_path):
        return None

    def extract_genome_name(path):
        """Extract genome name from BiG-SCAPE GBK path."""
        parts = path.split('/')
        for i, part in enumerate(parts):
            if part == 'antismash_input' and i + 1 < len(parts):
                return parts[i + 1]
        return Path(path).parent.name

    try:
        conn = sqlite3.connect(bigscape_db_path)
        cursor = conn.cursor()

        # Get genome -> GCF mapping
        cursor.execute("""
            SELECT g.path, bf.family_id
            FROM gbk g
            JOIN bgc_record b ON g.id = b.gbk_id
            LEFT JOIN bgc_record_family bf ON b.id = bf.record_id
            WHERE bf.family_id IS NOT NULL
        """)

        genome_gcfs = defaultdict(set)
        for path, family_id in cursor.fetchall():
            genome = extract_genome_name(path)
            genome_gcfs[genome].add(family_id)

        # Get BGC type counts for top types
        cursor.execute("""
            SELECT product, COUNT(*) as count
            FROM bgc_record
            GROUP BY product
            HAVING count >= 50
            ORDER BY count DESC
            LIMIT 6
        """)
        top_types = [(row[0], row[1]) for row in cursor.fetchall()]

        conn.close()

        genomes = list(genome_gcfs.keys())
        n_genomes = len(genomes)
        total_gcfs = len(set().union(*genome_gcfs.values())) if genome_gcfs else 0

        if n_genomes < 10:
            return None

        # Calculate rarefaction curve
        n_points = min(50, n_genomes)
        x_values = sorted(set(np.linspace(1, n_genomes, n_points).astype(int)))

        all_curves = []
        for _ in range(n_iterations):
            shuffled = genomes.copy()
            random.shuffle(shuffled)
            seen_gcfs = set()
            curve = []
            for i, genome in enumerate(shuffled, 1):
                seen_gcfs.update(genome_gcfs[genome])
                if i in x_values:
                    curve.append(len(seen_gcfs))
            all_curves.append(curve)

        all_curves = np.array(all_curves)
        mean_gcfs = np.mean(all_curves, axis=0)
        lower_ci = np.percentile(all_curves, 2.5, axis=0)
        upper_ci = np.percentile(all_curves, 97.5, axis=0)

        # Calculate saturation
        if len(mean_gcfs) >= 10:
            early_rate = (mean_gcfs[len(mean_gcfs)//10] - mean_gcfs[0]) / (x_values[len(x_values)//10] - x_values[0]) if x_values[len(x_values)//10] != x_values[0] else 0
            late_rate = (mean_gcfs[-1] - mean_gcfs[-len(mean_gcfs)//10]) / (x_values[-1] - x_values[-len(x_values)//10]) if x_values[-1] != x_values[-len(x_values)//10] else 0
            saturation = (1 - (late_rate / early_rate)) * 100 if early_rate > 0 else 100
        else:
            saturation = 0

        # Plot
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
        ax.set_facecolor('white')

        ax.plot(x_values, mean_gcfs, color='#2c5aa0', linewidth=2.5, label='All BGC types')
        ax.fill_between(x_values, lower_ci, upper_ci, color='#2c5aa0', alpha=0.2)

        ax.set_xlabel('Number of Genomes Sampled', fontsize=12)
        ax.set_ylabel('Unique Gene Cluster Families (GCFs)', fontsize=12)
        ax.set_title(f'GCF Rarefaction Curve - {taxon}', fontsize=14)
        ax.grid(True, alpha=0.3)

        # Add annotation
        ax.text(0.02, 0.98,
                f'Saturation: {saturation:.0f}%\n{total_gcfs} GCFs from {n_genomes} genomes',
                transform=ax.transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='#e8f4e8', alpha=0.8, edgecolor='#4a9'))

        # Add interpretation guide
        ax.text(0.98, 0.02,
                'Plateau = diversity saturated\nRising = more diversity to discover',
                transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
                ha='right', color='#666')

        plt.tight_layout()
        output_path = Path(outdir) / 'rarefaction_curve.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        return {
            'generated': True,
            'n_genomes': n_genomes,
            'total_gcfs': total_gcfs,
            'saturation': saturation,
            'top_types': top_types
        }

    except Exception as e:
        print(f"Warning: Could not generate rarefaction curve: {e}")
        return None


def generate_bigslice_stats_html(bigslice_stats_file):
    '''Generate HTML for BiG-SLiCE clustering statistics'''
    import os

    if not os.path.exists(bigslice_stats_file) or os.path.basename(bigslice_stats_file).startswith('NO_'):
        return ''

    try:
        with open(bigslice_stats_file, 'r') as f:
            stats = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read BiG-SLiCE statistics: {e}")
        return ''

    if 'error' in stats:
        return f'''
    <h3>BiG-SLiCE Clustering Statistics</h3>
    <div class="info-box" style="background-color: #fff3cd; border-left: 4px solid #ffc107;">
        <p><strong>Note:</strong> {stats['error']}</p>
    </div>
    '''

    # Generate fragmentation breakdown
    frag_html = ''
    if 'fragmentation' in stats:
        frag_rows = []
        for frag_type, frag_data in stats['fragmentation'].items():
            frag_rows.append(f'''
                <tr>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd;">{frag_type.capitalize()}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{frag_data['count']}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{frag_data['avg_distance']}</td>
                </tr>
            ''')
        frag_html = f'''
            <tr>
                <td colspan="3" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">Fragmentation Breakdown:</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold;">Type</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">Count</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">Avg Distance</td>
            </tr>
            {''.join(frag_rows)}
        '''

    # Generate GCF size distribution preview (top 5)
    gcf_dist_html = ''
    if 'gcf_distribution' in stats and stats['gcf_distribution']:
        top_gcfs = stats['gcf_distribution'][:5]
        gcf_rows = []
        for gcf in top_gcfs:
            gcf_rows.append(f'''
                <tr>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd;">GCF {gcf['gcf_id']}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{gcf['bgc_count']}</td>
                </tr>
            ''')
        more_text = ''
        if len(stats['gcf_distribution']) > 5:
            more_text = f'''
                <tr>
                    <td colspan="2" style="padding: 8px; text-align: center; color: #666; font-style: italic;">
                        ... and {len(stats['gcf_distribution']) - 5} more GCFs
                    </td>
                </tr>
            '''
        gcf_dist_html = f'''
            <tr>
                <td colspan="2" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">Top Gene Cluster Families:</td>
            </tr>
            {''.join(gcf_rows)}
            {more_text}
        '''

    return f'''
    <h3>BiG-SLiCE Clustering Statistics</h3>
    <div class="info-box">
        <table style="width: 100%; border-collapse: collapse;">
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; font-weight: bold;">Total BGCs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_bgcs', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; font-weight: bold;">Total GCFs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_gcfs', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Distance Threshold (T)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('threshold', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Avg. Distance to GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('avg_distance_to_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">BGCs Assigned (d ≤ T)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('bgcs_assigned', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">BGCs Not Assigned (d > T)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('bgcs_not_assigned', 0)}</td>
            </tr>
            <tr>
                <td colspan="2" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">GCF Size Statistics:</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Min BGCs per GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('min_bgcs_per_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Avg BGCs per GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('avg_bgcs_per_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Max BGCs per GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('max_bgcs_per_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Singleton GCFs (n=1)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('singleton_gcfs', 0)}</td>
            </tr>
            {frag_html}
            {gcf_dist_html}
        </table>
        <p style="margin-top: 15px; font-size: 0.9em; color: #666;">
            <strong>Note:</strong> Statistics extracted from BiG-SLiCE database. For detailed analysis, use the BiG-SLiCE web interface or query the database directly.
        </p>
    </div>
    '''

def generate_bigscape_stats_html(bigscape_stats_file, mibig_included=False):
    '''Generate HTML for BiG-SCAPE clustering statistics'''
    import os

    if not os.path.exists(bigscape_stats_file) or os.path.basename(bigscape_stats_file).startswith('NO_'):
        return ''

    try:
        with open(bigscape_stats_file, 'r') as f:
            stats = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read BiG-SCAPE statistics: {e}")
        return ''

    if 'error' in stats:
        return f'''
    <h3>BiG-SCAPE Clustering Statistics</h3>
    <div class="info-box" style="background-color: #fff3cd; border-left: 4px solid #ffc107;">
        <p><strong>Note:</strong> {stats['error']}</p>
    </div>
    '''

    # Generate BGC class breakdown if available
    class_breakdown_html = ''
    if 'bgc_classes' in stats and stats['bgc_classes']:
        class_rows = []
        for bgc_class, class_data in sorted(stats['bgc_classes'].items()):
            class_rows.append(f'''
                <tr>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd;">{bgc_class}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{class_data['families']}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{class_data['bgcs']}</td>
                </tr>
            ''')
        class_breakdown_html = f'''
            <tr>
                <td colspan="2" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">BGC Class Breakdown:</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold;">Class</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">Families</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">BGCs</td>
            </tr>
            {''.join(class_rows)}
        '''

    # Build the complete statistics table
    cutoff_display = stats.get('cutoff', 'N/A')

    # Show N/A for MIBiG families if MIBiG was not included in the analysis
    if mibig_included:
        mibig_families_display = stats.get('families_with_mibig', 0)
    else:
        mibig_families_display = '<span style="color: #999;" title="MIBiG references were not included in this analysis">N/A</span>'

    return f'''
    <h3>BiG-SCAPE Clustering Statistics</h3>
    <div class="info-box">
        <table style="width: 100%; border-collapse: collapse;">
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; font-weight: bold;">Network Overview</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;"></td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Number of Families</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_families', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Total BGCs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_bgcs', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Average BGCs per Family</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('avg_bgcs_per_family', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Max BGCs in a Family</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('max_bgcs_per_family', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Singleton Families (n=1)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('singleton_families', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Families with MIBiG Reference BGCs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{mibig_families_display}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Distance Cutoff</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{cutoff_display}</td>
            </tr>
            {class_breakdown_html}
        </table>
        <p style="margin-top: 15px; font-size: 0.9em; color: #666;">
            <strong>Note:</strong> Statistics extracted from BiG-SCAPE clustering files. For detailed network analysis, open the BiG-SCAPE HTML output in <code>results/bigscape_results/[taxon]/index.html</code>.
        </p>
    </div>
    '''


def extract_assembly_id_from_genome_name(genome_name):
    """Extract GCA/GCF assembly ID from genome name."""
    import re
    # Look for GCA_XXXXXXXXX.X or GCF_XXXXXXXXX.X pattern
    match = re.search(r'(GC[AF]_\d+\.\d+)', genome_name)
    if match:
        return match.group(1)
    return None


def build_gcf_taxonomy_distribution(gcf_data, taxonomy_map, name_map=None):
    """
    Build GCF x Taxonomy distribution data.

    Args:
        gcf_data: GCF data from extract_gcf_representatives.py
        taxonomy_map: Taxonomy map keyed by assembly ID
        name_map: Optional dict mapping assembly_id -> genome_name (used to create reverse lookup)

    Returns:
        dict with:
        - 'heatmap_data': {gcf_id: {genus: count, ...}, ...}
        - 'gcf_taxonomy': {gcf_id: {'genera': {genus: count}, 'species': {species: count}}, ...}
        - 'genera': sorted list of all genera
        - 'species': sorted list of all species
        - 'top_gcfs': list of top GCFs by member count
    """
    if not gcf_data or not taxonomy_map:
        return None

    genome_gcf_mapping = gcf_data.get('genome_gcf_mapping', {})
    family_metadata = gcf_data.get('family_metadata', {})
    gcfs = gcf_data.get('gcfs', [])

    if not genome_gcf_mapping:
        return None

    # Enhance family_metadata with KCB hit info from gcfs array
    for gcf in gcfs:
        family_id = str(gcf.get('family_id', ''))
        if family_id in family_metadata:
            kcb_hit = gcf.get('kcb_hit', '')
            kcb_acc = gcf.get('kcb_acc', '')
            # Handle NaN/float values
            if kcb_hit is None or (isinstance(kcb_hit, float) and str(kcb_hit) == 'nan'):
                kcb_hit = ''
            if kcb_acc is None or (isinstance(kcb_acc, float) and str(kcb_acc) == 'nan'):
                kcb_acc = ''
            family_metadata[family_id]['kcb_hit'] = str(kcb_hit) if kcb_hit else ''
            family_metadata[family_id]['kcb_acc'] = str(kcb_acc) if kcb_acc else ''

    # Build reverse name_map: genome_name -> assembly_id
    reverse_name_map = {}
    if name_map:
        reverse_name_map = {v: k for k, v in name_map.items()}

    # Build genome -> taxonomy lookup
    genome_to_taxonomy = {}
    for genome_name in genome_gcf_mapping.keys():
        # First try using reverse name_map
        assembly_id = reverse_name_map.get(genome_name)
        # Fallback: try extracting from genome name (for names that include assembly ID)
        if not assembly_id:
            assembly_id = extract_assembly_id_from_genome_name(genome_name)
        if assembly_id and assembly_id in taxonomy_map:
            tax_info = taxonomy_map[assembly_id]
            lineage = tax_info.get('lineage', {})
            genus = lineage.get('genus', {}).get('name', 'Unknown')
            species = lineage.get('species', {}).get('name', 'Unknown')
            genome_to_taxonomy[genome_name] = {
                'genus': genus,
                'species': species
            }

    # Build GCF -> taxonomy distribution
    gcf_taxonomy = {}  # {gcf_id: {'genera': {genus: count}, 'species': {species: count}}}
    all_genera = set()
    all_species = set()

    for genome_name, gcf_ids in genome_gcf_mapping.items():
        tax = genome_to_taxonomy.get(genome_name)
        if not tax:
            continue

        genus = tax['genus']
        species = tax['species']
        all_genera.add(genus)
        all_species.add(species)

        for gcf_id in gcf_ids:
            gcf_id_str = str(gcf_id)
            if gcf_id_str not in gcf_taxonomy:
                gcf_taxonomy[gcf_id_str] = {'genera': {}, 'species': {}}

            gcf_taxonomy[gcf_id_str]['genera'][genus] = gcf_taxonomy[gcf_id_str]['genera'].get(genus, 0) + 1
            gcf_taxonomy[gcf_id_str]['species'][species] = gcf_taxonomy[gcf_id_str]['species'].get(species, 0) + 1

    # Build heatmap data (GCF x Genus)
    heatmap_data = {}
    for gcf_id, tax_dist in gcf_taxonomy.items():
        heatmap_data[gcf_id] = tax_dist['genera']

    # Get top GCFs by member count
    top_gcfs = sorted(
        [(gcf_id, family_metadata.get(gcf_id, {}).get('member_count', 0), family_metadata.get(gcf_id, {}).get('product', 'Unknown'))
         for gcf_id in gcf_taxonomy.keys()],
        key=lambda x: -x[1]
    )[:50]  # Top 50 GCFs

    return {
        'heatmap_data': heatmap_data,
        'gcf_taxonomy': gcf_taxonomy,
        'genera': sorted(all_genera),
        'species': sorted(all_species),
        'top_gcfs': top_gcfs,
        'family_metadata': family_metadata
    }


def generate_bgc_distribution_html(gcf_data, taxonomy_map, gtdbtk_summary_path=None, name_map=None):
    """
    Generate HTML for BGC Distribution tab (replaces Tree View).
    Shows GCF x Taxonomy heatmap and distribution insights.
    """
    # Try to use GTDB-Tk taxonomy if available (more accurate)
    taxonomy_source = "NCBI"
    if gtdbtk_summary_path:
        try:
            gtdbtk_df = pd.read_csv(gtdbtk_summary_path, sep='\t')
            # Build taxonomy map from GTDB-Tk classification
            gtdb_taxonomy = {}
            for _, row in gtdbtk_df.iterrows():
                user_genome = row.get('user_genome', '')
                classification = row.get('classification', '')
                if user_genome and classification:
                    # Parse GTDB classification string: d__Bacteria;p__Proteobacteria;c__...;g__Genus;s__Species
                    parts = classification.split(';')
                    lineage = {}
                    for part in parts:
                        if part.startswith('g__'):
                            lineage['genus'] = part[3:] or 'Unknown'
                        elif part.startswith('s__'):
                            lineage['species'] = part[3:] or 'Unknown'

                    assembly_id = extract_assembly_id_from_genome_name(user_genome)
                    if assembly_id:
                        gtdb_taxonomy[assembly_id] = {
                            'lineage': {
                                'genus': {'name': lineage.get('genus', 'Unknown')},
                                'species': {'name': lineage.get('species', 'Unknown')}
                            }
                        }

            if gtdb_taxonomy:
                taxonomy_map = gtdb_taxonomy
                taxonomy_source = "GTDB-Tk"
                print(f"Using GTDB-Tk taxonomy for {len(gtdb_taxonomy)} genomes")
        except Exception as e:
            print(f"Warning: Could not load GTDB-Tk taxonomy: {e}")

    dist_data = build_gcf_taxonomy_distribution(gcf_data, taxonomy_map, name_map)

    if not dist_data:
        return '''
        <div class="info-box warning">
            <p><strong>BGC Distribution analysis requires:</strong></p>
            <ul>
                <li>BiG-SCAPE clustering results (GCF assignments)</li>
                <li>Taxonomy data (from NCBI or GTDB-Tk)</li>
            </ul>
            <p>Run with <code>--clustering bigscape</code> to enable GCF analysis.</p>
        </div>
        '''

    genera = dist_data['genera']
    top_gcfs = dist_data['top_gcfs']
    gcf_taxonomy = dist_data['gcf_taxonomy']
    family_metadata = dist_data['family_metadata']

    # Limit to top genera by total BGC count
    genus_totals = {}
    for gcf_id, tax_dist in gcf_taxonomy.items():
        for genus, count in tax_dist['genera'].items():
            genus_totals[genus] = genus_totals.get(genus, 0) + count

    top_genera = sorted(genus_totals.items(), key=lambda x: -x[1])[:20]
    top_genera_names = [g[0] for g in top_genera]

    # Build heatmap data for JavaScript
    heatmap_rows = []
    for gcf_id, member_count, product in top_gcfs[:30]:  # Top 30 GCFs for heatmap
        row_data = []
        for genus in top_genera_names:
            count = gcf_taxonomy.get(gcf_id, {}).get('genera', {}).get(genus, 0)
            row_data.append(count)
        heatmap_rows.append({
            'gcf_id': gcf_id,
            'product': product[:30],
            'member_count': member_count,
            'values': row_data
        })

    heatmap_json = json.dumps({
        'rows': heatmap_rows,
        'columns': top_genera_names
    })

    # Calculate distribution statistics
    gcf_specificity = []  # GCFs that are genus-specific
    gcf_widespread = []   # GCFs found in many genera

    for gcf_id, tax_dist in gcf_taxonomy.items():
        genera_with_gcf = [g for g, c in tax_dist['genera'].items() if c > 0]
        num_genera = len(genera_with_gcf)
        total_bgcs = sum(tax_dist['genera'].values())

        meta = family_metadata.get(gcf_id, {})
        product = meta.get('product', 'Unknown')
        kcb_hit = meta.get('kcb_hit', '')
        kcb_acc = meta.get('kcb_acc', '')

        if num_genera == 1 and total_bgcs >= 3:
            gcf_specificity.append({
                'gcf_id': gcf_id,
                'genus': genera_with_gcf[0],
                'count': total_bgcs,
                'product': product,
                'kcb_hit': kcb_hit,
                'kcb_acc': kcb_acc
            })
        elif num_genera >= 5:
            gcf_widespread.append({
                'gcf_id': gcf_id,
                'num_genera': num_genera,
                'count': total_bgcs,
                'product': product,
                'kcb_hit': kcb_hit,
                'kcb_acc': kcb_acc
            })

    # Sort by count
    gcf_specificity.sort(key=lambda x: -x['count'])
    gcf_widespread.sort(key=lambda x: -x['num_genera'])

    # Helper to create status badge HTML
    def make_status_badge(kcb_hit, kcb_acc):
        if kcb_hit:
            # Known compound - green badge with link to MIBiG if available
            hit_display = kcb_hit[:25] + '...' if len(kcb_hit) > 25 else kcb_hit
            if kcb_acc:
                return f'<a href="https://mibig.secondarymetabolites.org/repository/{kcb_acc}" target="_blank" style="background: #d4edda; color: #155724; padding: 2px 6px; border-radius: 3px; font-size: 0.85em; text-decoration: none;" title="{kcb_hit}">{hit_display}</a>'
            else:
                return f'<span style="background: #d4edda; color: #155724; padding: 2px 6px; border-radius: 3px; font-size: 0.85em;" title="{kcb_hit}">{hit_display}</span>'
        else:
            # Potentially novel - amber badge
            return '<span style="background: #fff3cd; color: #856404; padding: 2px 6px; border-radius: 3px; font-size: 0.85em;">Potentially Novel</span>'

    # Build specificity insights HTML
    specificity_rows = ''
    for item in gcf_specificity[:10]:
        status_badge = make_status_badge(item.get('kcb_hit', ''), item.get('kcb_acc', ''))
        specificity_rows += f'''
            <tr>
                <td style="padding: 6px 10px;">GCF-{item['gcf_id']}</td>
                <td style="padding: 6px 10px;">{item['product'][:30]}</td>
                <td style="padding: 6px 10px;">{status_badge}</td>
                <td style="padding: 6px 10px; font-weight: bold;">{item['genus']}</td>
                <td style="padding: 6px 10px; text-align: right;">{item['count']}</td>
            </tr>
        '''

    widespread_rows = ''
    for item in gcf_widespread[:10]:
        status_badge = make_status_badge(item.get('kcb_hit', ''), item.get('kcb_acc', ''))
        widespread_rows += f'''
            <tr>
                <td style="padding: 6px 10px;">GCF-{item['gcf_id']}</td>
                <td style="padding: 6px 10px;">{item['product'][:30]}</td>
                <td style="padding: 6px 10px;">{status_badge}</td>
                <td style="padding: 6px 10px; text-align: right;">{item['num_genera']}</td>
                <td style="padding: 6px 10px; text-align: right;">{item['count']}</td>
            </tr>
        '''

    return f'''
    <p style="color: #666; margin-bottom: 20px;">
        <em>Analysis of BGC distribution across taxonomic groups. Taxonomy source: <strong>{taxonomy_source}</strong>.
        Shows which Gene Cluster Families (GCFs) are taxon-specific vs widespread.</em>
    </p>

    <h3>GCF x Genus Heatmap</h3>
    <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
        Top 30 GCFs (rows) vs top 20 genera (columns). Color intensity = number of BGCs.
    </p>
    <div id="heatmap-container" style="width: 100%; max-width: 900px; margin-bottom: 30px;">
        <canvas id="heatmap-canvas"></canvas>
    </div>

    <script>
        (function() {{
            const data = {heatmap_json};
            const canvas = document.getElementById('heatmap-canvas');
            const ctx = canvas.getContext('2d');

            const rows = data.rows;
            const cols = data.columns;

            // Fixed dimensions for consistent sizing
            const labelWidth = 200;
            const rightMargin = 50;
            const cellWidth = 55;
            const cellHeight = 22;
            const labelHeight = 100;

            // Calculate canvas size (minimum width for legend text)
            const minWidth = 480;
            const width = Math.max(minWidth, labelWidth + cols.length * cellWidth + rightMargin);
            const height = labelHeight + rows.length * cellHeight + 40;

            // High-DPI canvas for sharp text
            const dpr = window.devicePixelRatio || 1;
            canvas.width = width * dpr;
            canvas.height = height * dpr;
            canvas.style.width = width + 'px';
            canvas.style.height = height + 'px';
            ctx.scale(dpr, dpr);

            // Find max value for color scaling
            let maxVal = 1;
            rows.forEach(row => {{
                row.values.forEach(v => {{ if (v > maxVal) maxVal = v; }});
            }});

            // Color scale function (white to blue)
            function getColor(value) {{
                if (value === 0) return '#f8f9fa';
                const intensity = Math.min(1, value / maxVal);
                const r = Math.round(255 - intensity * 212);
                const g = Math.round(255 - intensity * 165);
                const b = Math.round(255 - intensity * 95);
                return `rgb(${{r}}, ${{g}}, ${{b}})`;
            }}

            ctx.fillStyle = 'white';
            ctx.fillRect(0, 0, width, height);

            // Draw column labels (genera) - rotated
            ctx.save();
            ctx.font = 'bold 11px sans-serif';
            ctx.fillStyle = '#333';
            cols.forEach((col, i) => {{
                ctx.save();
                ctx.translate(labelWidth + i * cellWidth + cellWidth/2, labelHeight - 5);
                ctx.rotate(-Math.PI / 4);
                ctx.textAlign = 'left';
                ctx.fillText(col.length > 15 ? col.slice(0, 15) + '...' : col, 0, 0);
                ctx.restore();
            }});
            ctx.restore();

            // Draw rows
            rows.forEach((row, i) => {{
                const y = labelHeight + i * cellHeight;

                // Row label (GCF)
                ctx.font = '11px sans-serif';
                ctx.fillStyle = '#333';
                ctx.textAlign = 'right';
                const label = `GCF-${{row.gcf_id}} ${{row.product}}`;
                ctx.fillText(label.length > 26 ? label.slice(0, 26) + '...' : label, labelWidth - 5, y + cellHeight/2 + 4);

                // Cells
                row.values.forEach((value, j) => {{
                    const x = labelWidth + j * cellWidth;
                    ctx.fillStyle = getColor(value);
                    ctx.fillRect(x, y, cellWidth - 1, cellHeight - 1);

                    // Show value if > 0
                    if (value > 0) {{
                        ctx.fillStyle = value > maxVal * 0.5 ? 'white' : '#333';
                        ctx.font = 'bold 10px sans-serif';
                        ctx.textAlign = 'center';
                        ctx.fillText(value.toString(), x + cellWidth/2, y + cellHeight/2 + 3);
                    }}
                }});

                // Member count on right
                ctx.fillStyle = '#666';
                ctx.font = '10px sans-serif';
                ctx.textAlign = 'left';
                ctx.fillText(`(${{row.member_count}})`, labelWidth + cols.length * cellWidth + 5, y + cellHeight/2 + 3);
            }});

            // Legend
            const legendY = labelHeight + rows.length * cellHeight + 25;
            ctx.font = '11px sans-serif';
            ctx.fillStyle = '#666';
            ctx.textAlign = 'left';
            ctx.fillText('Color intensity: BGC count (darker = more BGCs in that genus)', 10, legendY);
        }})();
    </script>

    <div style="display: flex; gap: 30px; flex-wrap: wrap; margin-top: 20px;">
        <div style="flex: 1; min-width: 450px;">
            <h3>Genus-Specific GCFs</h3>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
                GCFs found in only one genus (potential taxon-specific metabolites)
            </p>
            <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                <thead>
                    <tr style="background: #f0f4f8;">
                        <th style="padding: 8px 10px; text-align: left;">GCF</th>
                        <th style="padding: 8px 10px; text-align: left;">Product</th>
                        <th style="padding: 8px 10px; text-align: left;">Status</th>
                        <th style="padding: 8px 10px; text-align: left;">Genus</th>
                        <th style="padding: 8px 10px; text-align: right;">BGCs</th>
                    </tr>
                </thead>
                <tbody>
                    {specificity_rows if specificity_rows else '<tr><td colspan="5" style="padding: 10px; color: #666;">No genus-specific GCFs found</td></tr>'}
                </tbody>
            </table>
        </div>

        <div style="flex: 1; min-width: 450px;">
            <h3>Widespread GCFs</h3>
            <p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">
                GCFs found across 5+ genera (conserved or horizontally transferred)
            </p>
            <table style="width: 100%; border-collapse: collapse; font-size: 0.9em;">
                <thead>
                    <tr style="background: #f0f4f8;">
                        <th style="padding: 8px 10px; text-align: left;">GCF</th>
                        <th style="padding: 8px 10px; text-align: left;">Product</th>
                        <th style="padding: 8px 10px; text-align: left;">Status</th>
                        <th style="padding: 8px 10px; text-align: right;">Genera</th>
                        <th style="padding: 8px 10px; text-align: right;">BGCs</th>
                    </tr>
                </thead>
                <tbody>
                    {widespread_rows if widespread_rows else '<tr><td colspan="5" style="padding: 10px; color: #666;">No widespread GCFs found</td></tr>'}
                </tbody>
            </table>
        </div>
    </div>

    <div style="margin-top: 20px; padding: 12px 16px; background: linear-gradient(135deg, rgba(106, 90, 160, 0.1) 0%, rgba(78, 115, 180, 0.1) 100%); border: 1px solid rgba(106, 90, 160, 0.2); border-radius: 6px;">
        <p style="margin: 0; color: #4a5568; font-size: 0.95em;">
            <strong style="color: #6a5aa0;">See the Clustering tab</strong> for detailed GCF analysis including gene diagrams, domain annotations, and links to antiSMASH results.
        </p>
    </div>

    <div style="margin-top: 20px; padding: 15px; background: #f8f9fa; border-radius: 6px;">
        <h4 style="margin-top: 0;">Summary Statistics</h4>
        <ul style="margin: 10px 0; color: #555;">
            <li><strong>{len(gcf_taxonomy)}</strong> GCFs analyzed across <strong>{len(genera)}</strong> genera</li>
            <li><strong>{len(gcf_specificity)}</strong> genus-specific GCFs (potential taxon markers)</li>
            <li><strong>{len(gcf_widespread)}</strong> widespread GCFs (found in 5+ genera)</li>
        </ul>
        <p style="color: #666; font-size: 0.9em; margin-bottom: 0;">
            <strong>Interpretation:</strong> Genus-specific GCFs may represent recently evolved or taxon-restricted biosynthetic capabilities.
            Widespread GCFs may indicate ancestral pathways or horizontal gene transfer events.
        </p>
    </div>
    '''


def generate_gcf_visualization_html(gcf_data_file, taxon):
    '''Generate HTML for GCF representative visualization in Clustering tab'''
    import os

    if not os.path.exists(gcf_data_file):
        return ''

    try:
        with open(gcf_data_file, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read GCF data: {e}")
        return ''

    if 'error' in data:
        return f'''
    <div class="gcf-section">
        <h4>Gene Cluster Family Representatives</h4>
        <div class="info-box" style="background-color: #fff3cd; border-left: 4px solid #ffc107;">
            <p><strong>Note:</strong> {data['error']}</p>
        </div>
    </div>
    '''

    gcfs = data.get('gcfs', [])
    summary = data.get('summary', {})

    if not gcfs:
        return '''
    <div class="gcf-section">
        <h4>Gene Cluster Family Representatives</h4>
        <p style="color: #666;">No gene cluster families found.</p>
    </div>
    '''

    total = summary.get('total', len(gcfs))
    singletons = summary.get('singletons', 0)
    clusters = summary.get('clusters', 0)

    # Build GCF cards HTML
    gcf_cards = ''
    for gcf in gcfs:
        family_id = gcf.get('family_id', '?')
        member_count = gcf.get('member_count', 0)
        is_singleton = gcf.get('is_singleton', False)
        product = gcf.get('product', 'Unknown')
        organism = gcf.get('organism', 'Unknown organism')
        genome_name = gcf.get('genome_name', 'unknown')
        region_number = gcf.get('region_number', 1)
        svg_diagram = gcf.get('svg_diagram', '')
        genes = gcf.get('genes', [])
        antismash_link = gcf.get('antismash_link', '')

        badge_class = 'singleton' if is_singleton else 'cluster'
        badge_text = 'Singleton' if is_singleton else f'{member_count} BGCs'

        # Novelty/KCB badge
        kcb_hit = gcf.get('kcb_hit', '')
        kcb_acc = gcf.get('kcb_acc', '')
        # Handle NaN/float values from pandas
        if kcb_hit is None or (isinstance(kcb_hit, float) and str(kcb_hit) == 'nan'):
            kcb_hit = ''
        else:
            kcb_hit = str(kcb_hit)
        if kcb_acc is None or (isinstance(kcb_acc, float) and str(kcb_acc) == 'nan'):
            kcb_acc = ''
        else:
            kcb_acc = str(kcb_acc)

        if kcb_hit:
            # Has KCB hit - show the known cluster name with MIBiG link
            mibig_link = f'https://mibig.secondarymetabolites.org/repository/{kcb_acc}' if kcb_acc else '#'
            kcb_badge = f'<a href="{mibig_link}" target="_blank" style="background: #27ae60; color: white; padding: 2px 8px; border-radius: 4px; font-size: 0.75em; text-decoration: none; margin-left: 8px;" title="KnownClusterBlast hit: {kcb_hit}">{kcb_hit[:25]}{"..." if len(kcb_hit) > 25 else ""}</a>'
        else:
            # No KCB hit - potentially novel
            kcb_badge = '<span style="background: #8e44ad; color: white; padding: 2px 8px; border-radius: 4px; font-size: 0.75em; margin-left: 8px;" title="No KnownClusterBlast hits - may represent a novel BGC">Potentially Novel</span>'

        # Check if domain hits are stored directly in genes (new approach)
        # or in enhanced_analysis (legacy approach)
        has_enhanced = False
        use_gene_level_hits = False

        # Check first gene for direct domain hits
        if genes and (genes[0].get('clusterhmmer_hits') is not None or genes[0].get('tigrfam_hits') is not None):
            use_gene_level_hits = True
            has_enhanced = any(
                gene.get('clusterhmmer_hits') or gene.get('tigrfam_hits')
                for gene in genes
            )

        # Fallback: build lookup dicts from enhanced_analysis (legacy)
        clusterhmmer_by_gene = {}
        tigrfam_by_gene = {}
        if not use_gene_level_hits:
            enhanced_analysis = gcf.get('enhanced_analysis')
            if enhanced_analysis:
                clusterhmmer_hits = enhanced_analysis.get('clusterhmmer_hits', [])
                tigrfam_hits = enhanced_analysis.get('tigrfam_hits', [])
                has_enhanced = bool(clusterhmmer_hits or tigrfam_hits)

                # Build lookup by gene for ClusterHmmer
                for hit in clusterhmmer_hits:
                    gene_id = hit.get('gene', '')
                    if gene_id:
                        if gene_id not in clusterhmmer_by_gene:
                            clusterhmmer_by_gene[gene_id] = []
                        clusterhmmer_by_gene[gene_id].append(hit)

                # Build lookup by gene for TIGRFam
                for hit in tigrfam_hits:
                    gene_id = hit.get('gene', '')
                    if gene_id:
                        if gene_id not in tigrfam_by_gene:
                            tigrfam_by_gene[gene_id] = []
                        tigrfam_by_gene[gene_id].append(hit)

        # Build gene table rows with integrated enhanced analysis
        gene_rows = ''
        for gene in genes[:30]:  # Limit to first 30 genes for display
            locus_tag = gene.get('locus_tag', '-')
            gene_name = gene.get('gene_name', '-') or '-'
            gene_product = gene.get('product', 'hypothetical protein')
            function = gene.get('function', 'other')
            color = gene.get('color', '#95a5a6')

            # Get ClusterHmmer domains for this gene
            if use_gene_level_hits:
                ch_hits = gene.get('clusterhmmer_hits', [])
            else:
                ch_hits = clusterhmmer_by_gene.get(locus_tag, [])
            # Get TIGRFam annotations for this gene
            if use_gene_level_hits:
                tf_hits = gene.get('tigrfam_hits', [])
            else:
                tf_hits = tigrfam_by_gene.get(locus_tag, [])

            # Enhance product annotation if it's generic
            # Priority: TIGRFam description > Pfam description > original
            generic_products = ['hypothetical protein', 'putative protein', 'unknown',
                               'uncharacterized protein', 'predicted protein', 'conserved protein',
                               'domain-containing protein', 'family protein']
            product_is_generic = any(gp in gene_product.lower() for gp in generic_products)
            enhanced_product = gene_product
            product_source = ''

            if product_is_generic or gene_product == '-':
                # Try TIGRFam first (more curated)
                if tf_hits:
                    best_tf = max(tf_hits, key=lambda x: x.get('score', 0))
                    tf_desc = best_tf.get('description', '')
                    if tf_desc and len(tf_desc) > 5:
                        # Clean up TIGRFam description (remove prefix like "TIGR02320: ")
                        if ': ' in tf_desc:
                            tf_desc = tf_desc.split(': ', 1)[-1]
                        enhanced_product = tf_desc
                        product_source = 'TF'
                # Fall back to Pfam if no TIGRFam
                if product_source == '' and ch_hits:
                    best_ch = max(ch_hits, key=lambda x: x.get('score', 0))
                    ch_desc = best_ch.get('description', '')
                    if ch_desc and len(ch_desc) > 5:
                        enhanced_product = ch_desc
                        product_source = 'Pfam'

            # Format product cell with source indicator
            if product_source:
                product_cell = f'{enhanced_product[:35]}{"..." if len(enhanced_product) > 35 else ""} <span style="background: #3498db; color: white; padding: 1px 4px; border-radius: 3px; font-size: 0.7em;" title="Inferred from {product_source}">{product_source}</span>'
            else:
                product_cell = f'{gene_product[:40]}{"..." if len(gene_product) > 40 else ""}'

            # Format ClusterHmmer column
            if ch_hits:
                ch_domains = [f"{h.get('hmm', '')} ({h.get('score', 0):.0f})" for h in ch_hits[:3]]
                clusterhmmer_col = '<br>'.join(ch_domains)
                if len(ch_hits) > 3:
                    clusterhmmer_col += f'<br><span style="color: #888;">+{len(ch_hits)-3} more</span>'
            else:
                clusterhmmer_col = '<span style="color: #ccc;">-</span>'

            # Format TIGRFam column
            if tf_hits:
                tf_annots = [f"{h.get('domain', '')} ({h.get('score', 0):.0f})" for h in tf_hits[:2]]
                tigrfam_col = '<br>'.join(tf_annots)
                if len(tf_hits) > 2:
                    tigrfam_col += f'<br><span style="color: #888;">+{len(tf_hits)-2} more</span>'
            else:
                tigrfam_col = '<span style="color: #ccc;">-</span>'

            # Format BlastP link
            translation = gene.get('translation', '')
            if translation:
                blastp_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&QUERY={urllib.parse.quote(translation)}&LINK_LOC=protein&PAGE_TYPE=BlastSearch"
                blastp_col = f'<a href="{blastp_url}" target="_blank" style="display: inline-block; background: #2c5aa0; color: white; padding: 2px 6px; border-radius: 3px; font-size: 0.75em; text-decoration: none;" title="Search NCBI BlastP">BlastP</a>'
            else:
                blastp_col = '<span style="color: #ccc;">-</span>'

            gene_rows += f'''
                <tr>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">
                        <span style="display: inline-block; width: 10px; height: 10px; background: {color}; border-radius: 2px; margin-right: 6px;"></span>
                        {locus_tag}
                    </td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">{gene_name}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">{product_cell}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">{function}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee; font-size: 0.8em;">{clusterhmmer_col}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee; font-size: 0.8em;">{tigrfam_col}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee; text-align: center;">{blastp_col}</td>
                </tr>
            '''

        more_genes_note = f'<p style="color: #666; font-size: 0.85em; margin-top: 5px;">Showing first 30 of {len(genes)} genes</p>' if len(genes) > 30 else ''

        # AntiSMASH link
        antismash_link_html = f'<a href="{antismash_link}" target="_blank" style="color: #2c5aa0;">View in antiSMASH</a>' if antismash_link else ''

        gcf_cards += f'''
        <div class="gcf-card" data-type="{badge_class}" data-size="{member_count}" data-product="{product}" data-novelty="{'novel' if not kcb_hit else 'known'}">
            <div class="gcf-header" onclick="toggleGCF('gcf_{family_id}')">
                <span class="gcf-title">GCF {family_id}: {product}</span>
                <span class="gcf-badge {badge_class}">{badge_text}</span>
                {kcb_badge}
                <span class="gcf-organism">{organism}</span>
                <span class="gcf-toggle" id="gcf_{family_id}_toggle">+</span>
            </div>
            <div class="gcf-content" id="gcf_{family_id}_content" style="display: none;">
                <div class="gene-diagram">
                    {svg_diagram if svg_diagram else '<p style="color: #999;">Gene diagram not available</p>'}
                </div>
                <div class="gene-legend">
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #e74c3c; border-radius: 2px;"></span> Core biosynthetic</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #e67e22; border-radius: 2px;"></span> Additional biosynthetic</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #27ae60; border-radius: 2px;"></span> Regulatory</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #3498db; border-radius: 2px;"></span> Transport</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #9b59b6; border-radius: 2px;"></span> Resistance</span>
                    <span><span style="display: inline-block; width: 12px; height: 12px; background: #95a5a6; border-radius: 2px;"></span> Other</span>
                </div>
                <div class="gene-table-container" style="margin-top: 15px; max-height: 400px; overflow-y: auto;">
                    <table class="gene-table" style="width: 100%; border-collapse: collapse; font-size: 0.85em;">
                        <thead>
                            <tr style="background: #f5f5f5;">
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Locus Tag</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Gene</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Product</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Function</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">ClusterHmmer</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">TIGRFam</th>
                                <th style="padding: 8px; text-align: center; border-bottom: 2px solid #ddd;">BlastP</th>
                            </tr>
                        </thead>
                        <tbody>
                            {gene_rows}
                        </tbody>
                    </table>
                    {more_genes_note}
                </div>
                <div class="gcf-links" style="margin-top: 15px; padding-top: 10px; border-top: 1px solid #eee;">
                    {antismash_link_html}
                </div>
            </div>
        </div>
        '''

    # Build complete GCF visualization HTML
    gcf_html = f'''
    <div class="gcf-section" style="margin-top: 30px;">
        <h4>Gene Cluster Family Representatives</h4>
        <p style="color: #666; margin-bottom: 15px;">
            Representative BGCs for each Gene Cluster Family identified by BiG-SCAPE. Click on a family to expand and view the gene diagram and annotations.
        </p>

        <div class="gcf-summary" style="display: flex; gap: 15px; margin-bottom: 20px;">
            <div style="background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center;">
                <div style="font-size: 1.5em; font-weight: bold; color: #2c5aa0;">{total}</div>
                <div style="font-size: 0.85em; color: #666;">Total GCFs</div>
            </div>
            <div style="background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center;">
                <div style="font-size: 1.5em; font-weight: bold; color: #27ae60;">{clusters}</div>
                <div style="font-size: 0.85em; color: #666;">Clusters (2+ BGCs)</div>
            </div>
            <div style="background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center;">
                <div style="font-size: 1.5em; font-weight: bold; color: #f39c12;">{singletons}</div>
                <div style="font-size: 0.85em; color: #666;">Singletons</div>
            </div>
        </div>

        <div class="gcf-controls" style="margin-bottom: 15px; display: flex; gap: 10px; align-items: center; flex-wrap: wrap;">
            <label style="font-size: 0.9em; color: #666;">Filter:</label>
            <select id="gcfFilter" onchange="filterGCFs()" style="padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px;">
                <option value="all">All GCFs</option>
                <option value="cluster">Clusters Only</option>
                <option value="singleton">Singletons Only</option>
            </select>
            <select id="gcfNoveltyFilter" onchange="filterGCFs()" style="padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px;">
                <option value="all">All (Novel + Known)</option>
                <option value="novel">Potentially Novel Only</option>
                <option value="known">Known (KCB Hits) Only</option>
            </select>
            <label style="font-size: 0.9em; color: #666; margin-left: 15px;">Sort:</label>
            <select id="gcfSort" onchange="sortGCFs(this.value)" style="padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px;">
                <option value="size">By Size (largest first)</option>
                <option value="product">By Product Type</option>
            </select>
            <button onclick="toggleAllGCFs()" style="margin-left: 15px; padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px; background: #f8f9fa; cursor: pointer;">Expand/Collapse All</button>
        </div>

        <div class="gcf-container" id="gcfContainer">
            {gcf_cards}
        </div>
    </div>

    <style>
        .gcf-card {{
            border: 1px solid #ddd;
            border-radius: 8px;
            margin-bottom: 10px;
            background: white;
            overflow: hidden;
        }}
        .gcf-header {{
            padding: 12px 15px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 12px;
            background: #fafafa;
            transition: background 0.2s;
        }}
        .gcf-header:hover {{
            background: #f0f0f0;
        }}
        .gcf-title {{
            font-weight: 600;
            color: #333;
            flex-grow: 1;
        }}
        .gcf-badge {{
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.8em;
            font-weight: 500;
        }}
        .gcf-badge.singleton {{
            background: #fff3cd;
            color: #856404;
        }}
        .gcf-badge.cluster {{
            background: #d4edda;
            color: #155724;
        }}
        .gcf-organism {{
            color: #666;
            font-size: 0.85em;
            font-style: italic;
        }}
        .gcf-toggle {{
            font-size: 1.2em;
            color: #999;
            font-weight: bold;
        }}
        .gcf-content {{
            padding: 15px;
            border-top: 1px solid #eee;
        }}
        .gene-diagram {{
            padding: 10px;
            background: #f8f9fa;
            border-radius: 4px;
            overflow-x: auto;
        }}
        .gene-legend {{
            margin-top: 10px;
            font-size: 0.8em;
            color: #666;
        }}
    </style>

    <script>
        function toggleGCF(gcfId) {{
            const content = document.getElementById(gcfId + '_content');
            const toggle = document.getElementById(gcfId + '_toggle');
            if (content.style.display === 'none') {{
                content.style.display = 'block';
                toggle.textContent = '-';
            }} else {{
                content.style.display = 'none';
                toggle.textContent = '+';
            }}
        }}

        function filterGCFs() {{
            const typeFilter = document.getElementById('gcfFilter').value;
            const noveltyFilter = document.getElementById('gcfNoveltyFilter').value;
            const cards = document.querySelectorAll('.gcf-card');
            cards.forEach(card => {{
                const matchesType = typeFilter === 'all' || card.dataset.type === typeFilter;
                const matchesNovelty = noveltyFilter === 'all' || card.dataset.novelty === noveltyFilter;
                card.style.display = (matchesType && matchesNovelty) ? 'block' : 'none';
            }});
        }}

        function sortGCFs(by) {{
            const container = document.getElementById('gcfContainer');
            const cards = Array.from(container.querySelectorAll('.gcf-card'));
            cards.sort((a, b) => {{
                if (by === 'size') {{
                    return parseInt(b.dataset.size) - parseInt(a.dataset.size);
                }}
                return a.dataset.product.localeCompare(b.dataset.product);
            }});
            cards.forEach(card => container.appendChild(card));
        }}

        function toggleAllGCFs() {{
            const contents = document.querySelectorAll('.gcf-content');
            const anyHidden = Array.from(contents).some(c => c.style.display === 'none');
            contents.forEach(content => {{
                content.style.display = anyHidden ? 'block' : 'none';
                const id = content.id.replace('_content', '_toggle');
                const toggle = document.getElementById(id);
                if (toggle) toggle.textContent = anyHidden ? '-' : '+';
            }});
        }}
    </script>
    '''

    return gcf_html


def generate_taxonomy_tree_html(taxonomy_tree_data):
    '''Generate interactive HTML for the taxonomy tree with genome-level data'''
    tree = taxonomy_tree_data.get('tree', {})
    metadata = taxonomy_tree_data.get('metadata', {})

    def render_genome_list(genomes, level):
        '''Render genome list under a species node with heatmap coloring'''
        if not genomes:
            return ''

        # Get all unique BGC types from genomes
        bgc_type_columns = set()
        for genome_data in genomes:
            bgc_types = genome_data.get('bgc_types', {})
            bgc_type_columns.update(bgc_types.keys())
        bgc_type_columns = sorted(bgc_type_columns)

        # Calculate global max count across ALL BGC types for consistent heatmap coloring
        global_bgc_max = 0
        for genome_data in genomes:
            bgc_types = genome_data.get('bgc_types', {})
            for count in bgc_types.values():
                if count > global_bgc_max:
                    global_bgc_max = count
        global_bgc_max = global_bgc_max if global_bgc_max > 0 else 1

        # Calculate max total BGCs for total column heatmap
        total_max = max([g.get('total_bgcs', 0) for g in genomes])
        total_max = total_max if total_max > 0 else 1

        def get_green_bg_color(count, max_count):
            '''Convert count to sequential green background color for Total BGCs'''
            if count == 0:
                return ''
            intensity = count / max_count
            if intensity <= 0.2:
                return '#e8f5e9'
            elif intensity <= 0.4:
                return '#a5d6a7'
            elif intensity <= 0.6:
                return '#66bb6a'
            elif intensity <= 0.8:
                return '#43a047'
            else:
                return '#2e7d32'

        def get_red_font_color(count, max_count):
            '''Convert count to sequential red font color for BGC types'''
            if count == 0:
                return ''
            intensity = count / max_count
            if intensity <= 0.2:
                return '#ffcdd2'
            elif intensity <= 0.4:
                return '#ef5350'
            elif intensity <= 0.6:
                return '#e53935'
            elif intensity <= 0.8:
                return '#c62828'
            else:
                return '#b71c1c'

        genome_list_html = '<div class="genome-list">'
        genome_list_html += '<table class="genome-table">'
        genome_list_html += '<tr><th>Genome</th><th>Total BGCs</th>'

        for bgc_type in bgc_type_columns:
            genome_list_html += f'<th>{bgc_type}</th>'
        genome_list_html += '</tr>'

        sorted_genomes = sorted(genomes, key=lambda x: x.get('name', ''))

        for genome_data in sorted_genomes:
            genome_name = genome_data.get('name', 'Unknown')
            total = genome_data.get('total_bgcs', 0)
            bgc_types = genome_data.get('bgc_types', {})

            genome_list_html += f'<tr><td><a href="genomes/{genome_name}.html" class="genome-link">{genome_name}</a></td>'

            total_bg_color = get_green_bg_color(total, total_max)
            total_style = f' style="background-color: {total_bg_color}; color: white; font-weight: bold;"' if total_bg_color else ' style="font-weight: bold;"'
            genome_list_html += f'<td class="count-cell"{total_style}>{total}</td>'

            for bgc_type in bgc_type_columns:
                count = bgc_types.get(bgc_type, 0)
                font_color = get_red_font_color(count, global_bgc_max)
                style = f' style="color: {font_color}; font-weight: bold;"' if font_color else ''
                genome_list_html += f'<td class="count-cell"{style}>{count if count > 0 else ""}</td>'

            genome_list_html += '</tr>'

        genome_list_html += '</table></div>'
        return genome_list_html

    def render_node(node, level=0):
        '''Recursively render tree nodes'''
        name = node.get('name', 'Unknown')
        rank = node.get('rank', '')
        stats = node.get('stats', {})
        children = node.get('children', {})
        genomes = node.get('genomes', [])

        genome_count = stats.get('genome_count', 0)
        total_bgcs = stats.get('total_bgcs', 0)
        avg_bgcs = stats.get('avg_bgcs_per_genome', 0)
        genomes_with_bgcs = stats.get('genomes_with_bgcs', 0)
        bgc_distribution = stats.get('bgc_type_distribution', {})

        bgc_types_html = ''
        if bgc_distribution:
            bgc_types_list = []
            for bgc_type, type_stats in list(bgc_distribution.items())[:3]:
                count = type_stats.get('total_count', 0)
                pct = type_stats.get('percentage_of_genomes', 0)
                bgc_types_list.append(f"{bgc_type}: {count} ({pct:.0f}%)")
            bgc_types_html = ', '.join(bgc_types_list)
            if len(bgc_distribution) > 3:
                bgc_types_html += f" (+{len(bgc_distribution)-3} more)"

        node_id = f"node_{abs(hash(name + rank + str(level)))}"
        has_children = len(children) > 0
        has_genomes = len(genomes) > 0 and rank == 'species'

        toggle_icon = '&#9660;' if (has_children or has_genomes) else '&#9679;'
        toggle_class = 'toggle' if (has_children or has_genomes) else 'leaf'

        html = f'''
        <div class="tree-node level-{level}">
            <div class="node-header {toggle_class}" onclick="toggleNode('{node_id}')">
                <span class="toggle-icon">{toggle_icon}</span>
                <span class="node-name"><strong>{name}</strong> <em>({rank})</em></span>
                <span class="node-stats">
                    {genome_count} genomes | {total_bgcs} BGCs | avg {avg_bgcs:.2f}
                    {(' | ' + bgc_types_html) if bgc_types_html else ''}
                </span>
            </div>'''

        if has_children or has_genomes:
            html += f'<div class="node-children" id="{node_id}">'

            for child_name, child_node in children.items():
                html += render_node(child_node, level + 1)

            if has_genomes:
                html += render_genome_list(genomes, level + 1)

            html += '</div>'

        html += '</div>'
        return html

    tree_html = f'''
    <div class="taxonomy-tree-container">
        <style>
            .taxonomy-tree-container {{
                background: white;
                padding: 15px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 15px 0;
                font-size: 0.85em;
            }}
            .tree-node {{
                margin-left: 12px;  /* Consistent indent since nodes are nested */
            }}
            .tree-node.level-0 {{ margin-left: 0px; }}  /* Root has no indent */
            .node-header {{
                padding: 4px 6px;
                margin: 1px 0;
                cursor: pointer;
                border-radius: 3px;
                transition: background 0.2s;
                line-height: 1.3;
            }}
            .node-header:hover {{
                background: #f0f0f0;
            }}
            .node-header.leaf {{
                cursor: default;
            }}
            .toggle-icon {{
                display: inline-block;
                width: 14px;
                font-size: 0.8em;
                color: #2c5aa0;
            }}
            .node-name {{
                margin-right: 10px;
            }}
            .node-stats {{
                color: #666;
                font-size: 0.85em;
            }}
            .node-children {{
                margin-top: 1px;
            }}
            .genome-list {{
                margin: 6px 0 6px 16px;
                padding: 8px;
                background: #f8f9fa;
                border-radius: 4px;
                border-left: 2px solid #5b8ac5;
                overflow-x: auto;
                max-width: calc(100% - 20px);
            }}
            .genome-table {{
                width: 100%;
                border-collapse: collapse;
                font-size: 0.8em;
                background: white;
            }}
            .genome-table th {{
                background: #2c5aa0;
                color: white;
                padding: 4px 6px;
                text-align: left;
                font-weight: bold;
                border: 1px solid #ddd;
            }}
            .genome-table td {{
                padding: 4px 6px;
                border: 1px solid #ddd;
                text-align: center;
                transition: opacity 0.2s;
            }}
            .genome-table td:first-child {{
                text-align: left;
                background: white !important;
            }}
            .genome-table tr:hover td {{
                opacity: 0.8;
            }}
            .genome-link {{
                color: #2c5aa0;
                text-decoration: none;
                font-weight: 500;
            }}
            .genome-link:hover {{
                color: #5b8ac5;
                text-decoration: underline;
            }}
            .count-cell {{
                text-align: center;
                color: #999;
            }}
        </style>
        <script>
            function toggleNode(nodeId) {{
                const element = document.getElementById(nodeId);
                if (element) {{
                    const header = element.previousElementSibling;
                    const icon = header.querySelector('.toggle-icon');
                    if (element.style.display === 'none') {{
                        element.style.display = 'block';
                        icon.innerHTML = '&#9660;';
                    }} else {{
                        element.style.display = 'none';
                        icon.innerHTML = '&#9654;';
                    }}
                }}
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                const nodeChildren = document.querySelectorAll('.node-children');
                nodeChildren.forEach((node, index) => {{
                    if (index > 0) {{
                        node.style.display = 'none';
                        const icon = node.previousElementSibling.querySelector('.toggle-icon');
                        if (icon) icon.innerHTML = '&#9654;';
                    }}
                }});
            }});
        </script>
        {render_node(tree)}
    </div>'''

    return tree_html

def generate_html_report(outdir, taxon, table_header, table_rows, stats, tree_html='',
                         bigslice_stats_html='', bigscape_stats_html='', gcf_visualization_html='',
                         donut_chart_generated=False, phylo_tree_generated=False, genome_table_html='',
                         resource_usage_html='', phylo_tree_data=None, gcf_data=None, taxonomy_map=None,
                         bigslice_section_html='', versions_data=None, rarefaction_stats=None,
                         gtdbtk_summary_path=None, name_map=None):
    '''Generate tab-based HTML report combining all visualizations'''

    # Clean taxon name for URLs - match Nextflow sanitizeTaxon function
    import re
    taxon_clean = re.sub(r'[^a-zA-Z0-9_]', '_', taxon)
    taxon_clean = re.sub(r'_+', '_', taxon_clean).strip('_')

    # Circular tree display - only show when GTDB-Tk phylogenetic tree is available
    tree_image = ''
    tree_title = ''
    tree_description = ''

    if phylo_tree_generated:
        tree_image = 'circular_tree.png'
        tree_title = 'Phylogenetic Tree (GTDB-Tk)'
        tree_description = '''Phylogenetic tree reconstructed using neighbor-joining from GTDB-Tk pairwise distances (120 marker genes).
        Branch lengths represent evolutionary distance. The tree is unrooted; displayed with arbitrary root for visualization.'''

    # Build KCB stats section for dashboard
    kcb_stats = stats.get('kcb_stats', {})
    kcb_section = ''
    kcb_mapping_section = ''
    novel_bgcs_tab_content = '''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to identify potentially novel BGCs.</p>'''
    if kcb_stats.get('total_regions', 0) > 0:
        sim_breakdown = kcb_stats.get('similarity_breakdown', {})
        sim_html = ', '.join([f"{k}: {v}" for k, v in sim_breakdown.items()]) if sim_breakdown else 'N/A'
        kcb_section = f'''
        <div class="stat-box">
            <div class="stat-value">{kcb_stats.get('hit_percentage', 0)}%</div>
            <div class="stat-label">BGCs with Known Cluster Hits</div>
            <div style="font-size: 0.8em; color: #666; margin-top: 5px;">({kcb_stats.get('regions_with_hits', 0)} of {kcb_stats.get('total_regions', 0)} regions)</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{kcb_stats.get('unique_known_clusters', 0)}</div>
            <div class="stat-label">Unique Known Clusters Matched</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{kcb_stats.get('novel_bgc_count', 0)}</div>
            <div class="stat-label">Potentially Novel BGCs</div>
        </div>'''

        # Build potentially novel BGCs summary (for Overview tab) and details (for separate tab)
        novel_bgcs = kcb_stats.get('novel_bgcs', [])
        novel_count = kcb_stats.get('novel_bgc_count', 0)

        # Build GCF lookup from gcf_data if available
        # Maps (genome, region_name) -> {family_id, member_count}
        gcf_lookup = {}
        has_gcf_data = False
        if gcf_data and isinstance(gcf_data, dict):
            # Use the pre-computed bgc_to_gcf mapping (keys are "genome|region_name" strings)
            bgc_to_gcf = gcf_data.get('bgc_to_gcf', {})
            if bgc_to_gcf:
                for key_str, info in bgc_to_gcf.items():
                    # Parse the key back to (genome, region_name)
                    parts = key_str.split('|')
                    if len(parts) == 2:
                        genome = parts[0]
                        region_name = parts[1]
                        gcf_lookup[(genome, region_name)] = info
                        has_gcf_data = True

        if novel_bgcs:
            # Group novel BGCs by product type for summary
            product_counts = {}
            for bgc in novel_bgcs:
                product = bgc.get('product', 'Unknown')
                product_counts[product] = product_counts.get(product, 0) + 1

            # Sort products by count
            sorted_products = sorted(product_counts.items(), key=lambda x: x[1], reverse=True)

            # Build summary rows for Overview tab
            summary_rows = ''
            for product, count in sorted_products[:20]:  # Show top 20 product types
                pct = round(count / novel_count * 100, 1) if novel_count > 0 else 0
                summary_rows += f'''
                <tr>
                    <td style="font-weight: bold;">{product}</td>
                    <td style="text-align: center;">{count}</td>
                    <td style="text-align: center;">{pct}%</td>
                </tr>'''

            # Build detailed table rows with antiSMASH links (for Novel BGCs tab)
            detail_rows = ''
            for bgc in novel_bgcs:
                genome = bgc.get('genome', 'unknown')
                region = bgc.get('region', '?')
                region_name = bgc.get('region_name', region)
                record_index = bgc.get('record_index', 1)
                product = bgc.get('product', 'Unknown')
                contig_edge = bgc.get('contig_edge', '')
                record_id = bgc.get('record_id', '')
                edge_badge = '<span style="background: #e74c3c; color: white; padding: 1px 5px; border-radius: 3px; font-size: 0.75em;">edge</span>' if str(contig_edge).lower() == 'true' else ''
                # Build antiSMASH region link - format is #r{record_index}c{region} (e.g., #r11c1 for region 1 on record 11)
                # Path: main_data_visualization/ -> main_analysis_results/taxon/ -> main_analysis_results/ -> results/ -> antismash_results/
                antismash_link = f'../../../antismash_results/{taxon_clean}/{genome}/index.html#r{record_index}c{region}'

                # Look up GCF info if available
                gcf_cell = ''
                if has_gcf_data:
                    # Use region_name (e.g., "40.1") for unique lookup
                    # Convert to string in case pandas read it as float
                    gcf_info = gcf_lookup.get((genome, str(region_name)), {})
                    if gcf_info:
                        fid = gcf_info.get('family_id', '')
                        mc = gcf_info.get('member_count', 1)
                        gcf_cell = f'<td style="text-align: center;"><span style="background: #3498db; color: white; padding: 2px 8px; border-radius: 4px; font-size: 0.85em;">GCF {fid}</span></td><td style="text-align: center;">{mc}</td>'
                    else:
                        gcf_cell = '<td style="text-align: center; color: #999;">-</td><td style="text-align: center; color: #999;">-</td>'

                detail_rows += f'''
                <tr>
                    <td><a href="genomes/{genome}.html" style="color: #2c5aa0;">{genome[:40]}{"..." if len(genome) > 40 else ""}</a></td>
                    <td style="text-align: center;"><a href="{antismash_link}" target="_blank" style="color: #28a745; font-weight: bold;">Region {region_name}</a></td>
                    <td>{product}</td>
                    <td style="text-align: center;">{edge_badge}</td>
                    {gcf_cell}
                </tr>'''

            # Novel BGCs summary removed from Overview - details available in Novel BGCs tab
            kcb_mapping_section = ''

            # Build table header - include GCF columns if clustering data available
            gcf_header = '<th>GCF Family</th><th>Members</th>' if has_gcf_data else ''
            gcf_description = ' When BiG-SCAPE clustering is enabled, the GCF (Gene Cluster Family) assignment shows how these novel BGCs group together.' if has_gcf_data else ''

            # Novel BGCs details section (for separate tab)
            novel_bgcs_tab_content = f'''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666; margin-bottom: 15px;">
                <em>These BGC regions did not return any hits from KnownClusterBlast (KCB) analysis against the MIBiG database,
                suggesting they may encode novel or uncharacterized biosynthetic pathways. Regions marked as "edge" are on contig
                boundaries and may be incomplete.{gcf_description}</em>
            </p>
            <div class="search-box">
                <input type="text" id="novelSearch" placeholder="Search novel BGCs..." onkeyup="filterNovelBGCs()">
            </div>
            <div class="table-container">
                <table id="novelTable">
                    <thead>
                        <tr>
                            <th>Genome</th>
                            <th>antiSMASH Region</th>
                            <th>Product Type</th>
                            <th>Contig Edge</th>
                            {gcf_header}
                        </tr>
                    </thead>
                    <tbody id="novelTableBody">
                        {detail_rows}
                    </tbody>
                </table>
            </div>'''
        elif kcb_stats.get('total_regions', 0) > 0:
            # All BGCs matched known clusters - no section needed on Overview
            kcb_mapping_section = ''
            novel_bgcs_tab_content = '''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666;">All detected BGC regions matched characterized clusters in the MIBiG database. No potentially novel BGCs identified.</p>'''
        else:
            novel_bgcs_tab_content = '''
            <h2>Potentially Novel BGCs</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to identify potentially novel BGCs.</p>'''

    # Build KCB Hits tab content (grouped by known cluster for scalability)
    kcb_hits_tab_content = '''
            <h2>KnownClusterBlast Hits</h2>
            <p style="color: #666;">No KnownClusterBlast data available. Run antiSMASH with <code>--antismash_cb_knownclusters true</code> to identify known cluster matches.</p>'''

    if kcb_stats.get('total_regions', 0) > 0:
        cluster_mapping = kcb_stats.get('cluster_mapping', [])
        regions_with_hits = kcb_stats.get('regions_with_hits', 0)
        unique_clusters = kcb_stats.get('unique_known_clusters', 0)
        sim_breakdown = kcb_stats.get('similarity_breakdown', {})

        # Build similarity breakdown display
        sim_badges = ''
        for sim_level, count in sim_breakdown.items():
            color = '#27ae60' if sim_level == 'high' else '#f39c12' if sim_level == 'medium' else '#95a5a6'
            sim_badges += f'<span style="background: {color}; color: white; padding: 4px 12px; border-radius: 4px; margin-right: 8px;">{sim_level}: {count}</span>'

        # Build known clusters table rows (grouped by cluster - scales with MIBiG size, not genome count)
        kcb_table_rows = ''
        for item in cluster_mapping:
            known_cluster = item['known_cluster']
            mibig_acc = item['mibig_acc']
            hit_count = item['count']
            regions = item['regions']

            # Get unique product types and genomes
            product_types = list(set(r['product'] for r in regions))
            products_display = ', '.join(product_types[:3])
            if len(product_types) > 3:
                products_display += f' (+{len(product_types) - 3} more)'

            # Build antiSMASH BGC region links
            bgc_links = []
            for r in regions[:5]:
                genome = r['genome']
                region = r['region']
                region_name = r.get('region_name', region)
                record_index = r.get('record_index', 1)
                antismash_link = f'../../../antismash_results/{taxon_clean}/{genome}/index.html#r{record_index}c{region}'
                bgc_links.append(f'<a href="{antismash_link}" target="_blank" style="color: #28a745;">{genome} Region {region_name}</a>')
            bgc_display = ', '.join(bgc_links)
            if len(regions) > 5:
                bgc_display += f' (+{len(regions) - 5} more)'

            # Similarity breakdown for this cluster
            sim_counts = {}
            for r in regions:
                sim = r.get('similarity', 'unknown')
                sim_counts[sim] = sim_counts.get(sim, 0) + 1
            sim_display = ' / '.join([f'{k}: {v}' for k, v in sim_counts.items()])

            # Determine row color based on highest similarity level
            if 'high' in sim_counts:
                row_bg = 'background: rgba(39, 174, 96, 0.15);'  # Light green
            elif 'medium' in sim_counts:
                row_bg = 'background: rgba(243, 156, 18, 0.15);'  # Light orange
            elif 'low' in sim_counts:
                row_bg = 'background: rgba(149, 165, 166, 0.15);'  # Light gray
            else:
                row_bg = ''

            # MIBiG link - use full accession with version and trailing slash
            mibig_link = f'https://mibig.secondarymetabolites.org/repository/{mibig_acc}/' if mibig_acc else '#'

            kcb_table_rows += f'''
                <tr style="{row_bg}">
                    <td><a href="{mibig_link}" target="_blank" style="color: #2c5aa0; font-weight: bold;">{known_cluster}</a></td>
                    <td><a href="{mibig_link}" target="_blank" style="color: #666;">{mibig_acc}</a></td>
                    <td style="text-align: center; font-weight: bold;">{hit_count}</td>
                    <td>{products_display}</td>
                    <td>{sim_display}</td>
                    <td>{bgc_display}</td>
                </tr>'''

        kcb_hits_tab_content = f'''
            <h2>KnownClusterBlast Hits</h2>
            <div style="margin-bottom: 20px;">
                <strong>Total Hits by Similarity:</strong> {sim_badges}
            </div>
            <div style="background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 6px; padding: 15px; margin-bottom: 20px;">
                <strong>Similarity Legend:</strong>
                <div style="margin-top: 10px; display: flex; gap: 20px; flex-wrap: wrap;">
                    <div><span style="background: #27ae60; color: white; padding: 2px 8px; border-radius: 3px;">high</span> &gt;75% sequence similarity to MIBiG reference</div>
                    <div><span style="background: #f39c12; color: white; padding: 2px 8px; border-radius: 3px;">medium</span> 50-75% sequence similarity</div>
                    <div><span style="background: #95a5a6; color: white; padding: 2px 8px; border-radius: 3px;">low</span> 15-50% sequence similarity</div>
                </div>
                <p style="margin: 10px 0 0 0; font-size: 0.9em; color: #666;">
                    The "Similarity" column shows how many regions matched each known cluster at different similarity levels (e.g., "high: 1 / low: 2" means 1 region matched with &gt;75% similarity and 2 regions matched with 15-50% similarity).
                </p>
            </div>
            <div class="search-box">
                <input type="text" id="kcbSearch" placeholder="Search known clusters..." onkeyup="filterKCBHits()">
            </div>
            <div class="table-container">
                <table id="kcbTable">
                    <thead>
                        <tr>
                            <th>Known Cluster</th>
                            <th>MIBiG ID</th>
                            <th>Hits</th>
                            <th>Product Types</th>
                            <th>Similarity</th>
                            <th>BGC Regions</th>
                        </tr>
                    </thead>
                    <tbody id="kcbTableBody">
                        {kcb_table_rows}
                    </tbody>
                </table>
            </div>
            <p style="color: #666; font-style: italic; margin-top: 15px;">
                Table shows {unique_clusters} unique known clusters from MIBiG database. Each row represents a characterized BGC that matched one or more regions in your dataset.
            </p>'''

    # =========================================================================
    # BUILD NEW OVERVIEW SECTIONS
    # =========================================================================

    # --- Quick Stats (merged into Summary Statistics in template) ---
    kcb_stats_data = stats.get('kcb_stats', {})
    quick_stats_html = ''  # Merged into Summary Statistics

    # --- Clustering Summary ---
    clustering_summary_html = ''
    if gcf_data and isinstance(gcf_data, dict) and gcf_data.get('gcfs'):
        gcfs = gcf_data['gcfs']
        gcf_summary = gcf_data.get('summary', {})
        total_gcfs = gcf_summary.get('total', len(gcfs))
        singletons = gcf_summary.get('singletons', sum(1 for g in gcfs if g.get('is_singleton')))
        clusters = total_gcfs - singletons
        largest_gcf = max(gcfs, key=lambda x: x.get('member_count', 0)) if gcfs else None

        clustering_summary_html = f'''
    <div class="clustering-summary-section" style="margin-top: 20px;">
        <h3 style="margin-bottom: 12px;">Clustering Summary (BiG-SCAPE)</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px;">
            <div style="background: rgba(106, 90, 160, 0.15); border: 1px solid rgba(106, 90, 160, 0.3); padding: 12px 15px; border-radius: 8px;">
                <div style="font-size: 1.6em; font-weight: bold; color: #6a5aa0;">{total_gcfs}</div>
                <div style="color: #555; font-size: 0.9em;">Gene Cluster Families</div>
            </div>
            <div style="background: rgba(128, 90, 160, 0.15); border: 1px solid rgba(128, 90, 160, 0.3); padding: 12px 15px; border-radius: 8px;">
                <div style="font-size: 1.6em; font-weight: bold; color: #805aa0;">{clusters}</div>
                <div style="color: #555; font-size: 0.9em;">Multi-member GCFs</div>
            </div>
            <div style="background: rgba(140, 110, 160, 0.15); border: 1px solid rgba(140, 110, 160, 0.3); padding: 12px 15px; border-radius: 8px;">
                <div style="font-size: 1.6em; font-weight: bold; color: #8c6ea0;">{singletons}</div>
                <div style="color: #555; font-size: 0.9em;">Singletons</div>
            </div>
            <div style="background: rgba(150, 120, 170, 0.15); border: 1px solid rgba(150, 120, 170, 0.3); padding: 12px 15px; border-radius: 8px;">
                <div style="font-size: 1.6em; font-weight: bold; color: #9678aa;">{largest_gcf.get("member_count", 0) if largest_gcf else 0}</div>
                <div style="color: #555; font-size: 0.9em;">Largest GCF Size</div>
                <div style="font-size: 0.8em; margin-top: 3px; color: #777;">{largest_gcf.get("product", "")[:25] if largest_gcf else ""}</div>
            </div>
        </div>
    </div>'''

    # --- Taxonomic Coverage ---
    taxonomic_coverage_html = ''
    if taxonomy_map and isinstance(taxonomy_map, dict):
        # Count unique values at each taxonomic level
        # Handle nested lineage structure: {genome_id: {lineage: {level: {name: "..."}}}}
        tax_levels = {}
        for genome_id, tax_info in taxonomy_map.items():
            if isinstance(tax_info, dict):
                lineage = tax_info.get('lineage', {})
                for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
                    level_data = lineage.get(level, {})
                    level_name = level_data.get('name', '') if isinstance(level_data, dict) else level_data
                    if level_name:
                        if level not in tax_levels:
                            tax_levels[level] = set()
                        tax_levels[level].add(level_name)

        if tax_levels:
            # Proper plural forms for taxonomic levels
            plurals = {
                'phylum': 'Phyla', 'class': 'Classes', 'order': 'Orders',
                'family': 'Families', 'genus': 'Genera', 'species': 'Species'
            }
            tax_stats = []
            for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
                if level in tax_levels:
                    count = len(tax_levels[level])
                    plural_label = plurals.get(level, level.title() + 's')
                    tax_stats.append(f'<div style="text-align: center;"><div style="font-size: 1.5em; font-weight: bold; color: #2c5aa0;">{count}</div><div style="color: #666;">{plural_label}</div></div>')

            taxonomic_coverage_html = f'''
    <div class="taxonomic-coverage-section" style="margin-top: 30px; background: #f0f7ff; padding: 20px; border-radius: 10px;">
        <h3>Taxonomic Coverage</h3>
        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(100px, 1fr)); gap: 15px; margin-top: 15px;">
            {''.join(tax_stats)}
        </div>
    </div>'''

    # Build BiG-SCAPE section - only show if BiG-SCAPE data is available
    bigscape_section_html = ''
    if bigscape_stats_html or gcf_visualization_html:
        bigscape_section_html = f'''
            <div class="clustering-section">
                <h3>BiG-SCAPE Gene Cluster Families</h3>
                <div class="info-box">
                    <p>BiG-SCAPE clusters biosynthetic gene clusters into gene cluster families based on sequence similarity.</p>
                    <p style="margin-top: 15px;"><strong>To view interactive results:</strong></p>
                    <ol style="margin: 10px 0 10px 20px; line-height: 1.8;">
                        <li>Open: <code>results/bigscape_results/{taxon_clean}/index.html</code></li>
                        <li>Select database: <code>results/bigscape_results/{taxon_clean}/{taxon_clean}.db</code></li>
                    </ol>
                </div>
                {bigscape_stats_html if bigscape_stats_html else ''}
                {gcf_visualization_html if gcf_visualization_html else ''}
            </div>'''

    # Build software versions HTML section
    versions_html = ''
    if versions_data:
        version_rows = []
        # Define display names for tools
        tool_names = {
            'antismash': 'antiSMASH',
            'bigscape': 'BiG-SCAPE',
            'bigslice': 'BiG-SLiCE',
            'gtdbtk': 'GTDB-Tk',
            'taxonkit': 'TaxonKit',
            'nextflow': 'Nextflow',
            'pipeline_version': 'Pipeline Version'
        }
        for tool, version in versions_data.items():
            display_name = tool_names.get(tool, tool.title())
            version_rows.append(f'''
                <tr>
                    <td style="padding: 8px 12px; border-bottom: 1px solid #eee; font-weight: 500;">{display_name}</td>
                    <td style="padding: 8px 12px; border-bottom: 1px solid #eee; font-family: monospace;">{version}</td>
                </tr>''')

        versions_html = f'''
            <div class="versions-section" style="margin-top: 30px;">
                <h3>Software Versions</h3>
                <p style="color: #666; margin-bottom: 15px;">
                    <em>Software versions used in this pipeline run for reproducibility.</em>
                </p>
                <table style="width: 100%; max-width: 500px; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">
                    <thead>
                        <tr style="background: #2c5aa0; color: white;">
                            <th style="padding: 10px 12px; text-align: left;">Software</th>
                            <th style="padding: 10px 12px; text-align: left;">Version</th>
                        </tr>
                    </thead>
                    <tbody>
                        {''.join(version_rows)}
                    </tbody>
                </table>
            </div>'''

    # Build rarefaction curve section for Overview tab
    rarefaction_section = ''
    if rarefaction_stats and rarefaction_stats.get('generated'):
        saturation = rarefaction_stats.get('saturation', 0)
        n_genomes = rarefaction_stats.get('n_genomes', 0)
        total_gcfs = rarefaction_stats.get('total_gcfs', 0)

        # Interpretation based on saturation level
        if saturation >= 90:
            interpretation = "GCF diversity is well-saturated. Additional sampling within this taxon is unlikely to discover many new gene cluster families."
            interpretation_color = "#27ae60"
        elif saturation >= 70:
            interpretation = "GCF diversity is moderately saturated. Some new gene cluster families may still be discovered with additional sampling."
            interpretation_color = "#f39c12"
        else:
            interpretation = "GCF diversity is not saturated. Significant new diversity may be discovered with additional sampling."
            interpretation_color = "#e74c3c"

        rarefaction_section = f'''
            <!-- GCF Rarefaction Curve -->
            <div class="plot" style="max-width: 800px; margin: 30px auto;">
                <h3 style="margin-top: 0;">GCF Rarefaction Curve</h3>
                <p style="color: #666; font-size: 0.9em; margin-bottom: 15px;">
                    Shows how the number of unique Gene Cluster Families (GCFs) increases as more genomes are sampled.
                    A plateau indicates sampling saturation - additional genomes would yield diminishing discovery of new GCFs.
                </p>
                <img src="rarefaction_curve.png" alt="GCF Rarefaction Curve" style="max-width: 100%; height: auto;">
                <div style="margin-top: 15px; padding: 15px; background: #f8f9fa; border-radius: 8px; border-left: 4px solid {interpretation_color};">
                    <strong style="color: {interpretation_color};">Saturation: {saturation:.0f}%</strong>
                    <span style="color: #666; margin-left: 15px;">({total_gcfs} GCFs from {n_genomes} genomes)</span>
                    <p style="margin: 10px 0 0 0; color: #555; font-size: 0.95em;">{interpretation}</p>
                </div>
            </div>'''

    # Build phylogenetic tree download info section
    phylo_tree_download_section = ''
    if phylo_tree_generated:
        phylo_tree_download_section = f'''
            <!-- Phylogenetic Tree Download Info -->
            <div style="margin-top: 30px; padding: 20px; background: linear-gradient(135deg, #f0f7ff 0%, #e8f4f8 100%); border: 1px solid #c5dbe8; border-radius: 10px;">
                <h3 style="margin-top: 0; color: #2c5aa0; display: flex; align-items: center; gap: 10px;">
                    <span style="font-size: 1.3em;">&#128193;</span> Phylogenetic Tree Files
                </h3>
                <p style="color: #555; margin-bottom: 15px;">
                    GTDB-Tk phylogenetic placement results are available for use with external visualization software (e.g., iTOL, FigTree, ggtree).
                </p>
                <div style="background: white; padding: 15px; border-radius: 6px; border: 1px solid #ddd;">
                    <div style="margin-bottom: 10px;">
                        <strong style="color: #333;">Pruned Newick tree:</strong>
                        <code style="background: #f5f5f5; padding: 3px 8px; border-radius: 4px; font-size: 0.9em; margin-left: 8px;">pruned_phylo_tree.nwk</code>
                    </div>
                    <p style="color: #666; font-size: 0.9em; margin: 0;">
                        This tree contains only the genomes from your analysis, pruned from the GTDB reference tree.
                        It can be loaded directly into tree visualization tools for further annotation and analysis.
                    </p>
                </div>
                <p style="color: #777; font-size: 0.85em; margin-top: 12px; margin-bottom: 0;">
                    <strong>Full results location:</strong> <code style="background: #f5f5f5; padding: 2px 6px; border-radius: 3px;">results/gtdbtk_results/{taxon_clean}/</code>
                </p>
            </div>'''

    # Build BGC Distribution section HTML (replaces Tree View)
    distribution_section = generate_bgc_distribution_html(gcf_data, taxonomy_map, gtdbtk_summary_path, name_map)

    # Use distribution_section for the "BGC Distribution" tab (replaces old tree view)
    tree_section = distribution_section  # Keep variable name for template compatibility

    # Build genome table section (avoid nested f-strings)
    genome_table_rows = genome_table_html if genome_table_html else '<tr><td colspan="6">No genome data available</td></tr>'

    html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>BGC Analysis Report - {taxon}</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px 40px; background: #f8f9fa; color: #333; }}
        h1 {{ color: #333; text-align: center; margin-bottom: 5px; }}
        h2 {{ color: #333; border-bottom: 2px solid #5b8ac5; padding-bottom: 10px; margin-top: 30px; }}
        h3 {{ color: #444; margin-top: 25px; }}
        .subtitle {{ text-align: center; color: #666; font-size: 1.1em; margin-bottom: 20px; }}

        /* Tab Styles */
        .tabs {{
            margin-top: 20px;
        }}
        .tabs input[type="radio"] {{
            display: none;
        }}
        .tabs label {{
            display: inline-block;
            padding: 12px 24px;
            background: #e9ecef;
            color: #666;
            cursor: pointer;
            border-radius: 8px 8px 0 0;
            margin-right: 4px;
            font-weight: 500;
            transition: all 0.2s;
            border: 1px solid #ddd;
            border-bottom: none;
        }}
        .tabs label:hover {{
            background: #dee2e6;
            color: #333;
        }}
        .tabs input[type="radio"]:checked + label {{
            background: white;
            color: #2c5aa0;
            border-color: #5b8ac5;
            font-weight: 600;
        }}
        .tab-content {{
            display: none;
            background: white;
            padding: 30px;
            border: 1px solid #ddd;
            border-radius: 0 8px 8px 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            min-height: 400px;
        }}
        #tab1:checked ~ #content1,
        #tab2:checked ~ #content2,
        #tab3:checked ~ #content3,
        #tab4:checked ~ #content4,
        #tab5:checked ~ #content5,
        #tab6:checked ~ #content6,
        #tab7:checked ~ #content7,
        #tab8:checked ~ #content8 {{
            display: block;
        }}

        /* Stats Dashboard */
        .stats-container {{
            display: grid;
            grid-template-columns: repeat(8, 1fr);
            gap: 8px;
            margin: 15px 0;
        }}
        @media (max-width: 1200px) {{
            .stats-container {{
                grid-template-columns: repeat(4, 1fr);
            }}
        }}
        @media (max-width: 768px) {{
            .stats-container {{
                grid-template-columns: repeat(2, 1fr);
            }}
        }}
        .stat-box {{
            background: #f8f9fa;
            padding: 12px 8px;
            border-radius: 6px;
            border: 1px solid #e9ecef;
            text-align: center;
        }}
        .stat-box.highlight {{
            /* Same as regular stat-box */
        }}
        .stat-value {{
            font-size: 1.4em;
            font-weight: bold;
            color: #2c5aa0;
        }}
        .stat-label {{
            color: #666;
            margin-top: 3px;
            font-size: 0.8em;
        }}

        /* Plots */
        .plot {{
            margin: 20px 0;
            text-align: center;
            background: #fafafa;
            padding: 20px;
            border-radius: 8px;
            border: 1px solid #eee;
        }}
        .plot img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
        }}
        .plot-row {{
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            justify-content: center;
        }}
        .plot-row .plot {{
            flex: 1;
            min-width: 300px;
            max-width: 500px;
        }}

        /* Tables */
        table {{
            border-collapse: collapse;
            width: 100%;
            background: white;
            font-size: 0.9em;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 10px 12px;
            text-align: left;
        }}
        th {{
            background: #2c5aa0;
            color: white;
            font-weight: 600;
            position: sticky;
            top: 0;
        }}
        tr:nth-child(even) {{ background-color: #f8f9fa; }}
        tr:hover {{ background-color: #e3f2fd; }}
        .table-container {{
            max-height: 500px;
            overflow-y: auto;
            border: 1px solid #ddd;
            border-radius: 8px;
        }}

        /* Search Box */
        .search-box {{
            margin-bottom: 15px;
        }}
        .search-box input {{
            width: 100%;
            max-width: 400px;
            padding: 10px 15px;
            font-size: 1em;
            border: 2px solid #ddd;
            border-radius: 6px;
            outline: none;
            transition: border-color 0.2s;
        }}
        .search-box input:focus {{
            border-color: #5b8ac5;
        }}

        /* Info Box */
        .info-box {{
            background: #e8f4f8;
            border-left: 4px solid #2c5aa0;
            padding: 20px;
            margin: 20px 0;
            border-radius: 0 8px 8px 0;
        }}
        .info-box.warning {{
            background: #fff8e1;
            border-left-color: #ffc107;
        }}

        /* Links & Buttons */
        a {{ color: #2c5aa0; text-decoration: none; }}
        a:hover {{ color: #5b8ac5; text-decoration: underline; }}
        .btn {{
            display: inline-block;
            padding: 10px 20px;
            background: #2c5aa0;
            color: white;
            border-radius: 6px;
            font-weight: 500;
            text-decoration: none;
            transition: background 0.2s;
        }}
        .btn:hover {{ background: #5b8ac5; color: white; text-decoration: none; }}
        code {{
            background: #f4f4f4;
            padding: 3px 8px;
            border-radius: 4px;
            font-family: 'Consolas', monospace;
            font-size: 0.9em;
            border: 1px solid #e0e0e0;
        }}

        /* Clustering sections */
        .clustering-section {{
            margin-bottom: 30px;
            padding-bottom: 30px;
            border-bottom: 1px solid #eee;
        }}
        .clustering-section:last-child {{
            border-bottom: none;
        }}
    </style>
</head>
<body>
    <h1>BGC Analysis Report</h1>
    <p class="subtitle">Taxon: <strong>{taxon}</strong></p>

    <div class="tabs">
        <input type="radio" id="tab1" name="tabs" checked>
        <label for="tab1">Overview</label>

        <input type="radio" id="tab2" name="tabs">
        <label for="tab2">Taxonomy</label>

        <input type="radio" id="tab3" name="tabs">
        <label for="tab3">BGC Distribution</label>

        <input type="radio" id="tab4" name="tabs">
        <label for="tab4">Genomes</label>

        <input type="radio" id="tab5" name="tabs">
        <label for="tab5">Clustering</label>

        <input type="radio" id="tab6" name="tabs">
        <label for="tab6">Novel BGCs</label>

        <input type="radio" id="tab7" name="tabs">
        <label for="tab7">KCB Hits</label>

        <input type="radio" id="tab8" name="tabs">
        <label for="tab8">Resources</label>

        <!-- Tab 1: Overview -->
        <div class="tab-content" id="content1">
            <h2>Summary Statistics</h2>
            <!-- Genome Statistics -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin-top: 12px;">
                <div style="background: rgba(44, 90, 160, 0.15); border: 1px solid rgba(44, 90, 160, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #2c5aa0;">{stats.get('total_genomes', 0)}</div>
                    <div style="color: #555; font-size: 0.9em;">Total Genomes</div>
                </div>
                <div style="background: rgba(58, 109, 153, 0.15); border: 1px solid rgba(58, 109, 153, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #3a6d99;">{stats.get('genomes_with_bgcs', 0)}</div>
                    <div style="color: #555; font-size: 0.9em;">Genomes with BGCs</div>
                </div>
                <div style="background: rgba(74, 122, 143, 0.15); border: 1px solid rgba(74, 122, 143, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #4a7a8f;">{stats.get('genomes_with_no_bgcs', 0)}</div>
                    <div style="color: #555; font-size: 0.9em;">Genomes without BGCs</div>
                </div>
            </div>
            <!-- BGC Statistics -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin-top: 12px;">
                <div style="background: rgba(61, 139, 139, 0.15); border: 1px solid rgba(61, 139, 139, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #3d8b8b;">{stats.get('total_bgcs', 0)}</div>
                    <div style="color: #555; font-size: 0.9em;">Total BGCs Found</div>
                </div>
                <div style="background: rgba(74, 149, 144, 0.15); border: 1px solid rgba(74, 149, 144, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #4a9590;">{stats.get('avg_bgcs', '0')}</div>
                    <div style="color: #555; font-size: 0.9em;">Avg BGCs/Genome</div>
                    <div style="font-size: 0.8em; margin-top: 3px; color: #777;">Range: {stats.get('min_bgcs', 0)} - {stats.get('max_bgcs', 0)}</div>
                </div>
                <div style="background: rgba(90, 159, 149, 0.15); border: 1px solid rgba(90, 159, 149, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #5a9f95;">{stats.get('most_common_bgc', 'N/A')}</div>
                    <div style="color: #555; font-size: 0.9em;">Most Common BGC Type</div>
                    <div style="font-size: 0.8em; margin-top: 3px; color: #777;">({stats.get('most_common_count', 0)} occurrences)</div>
                </div>
            </div>
            <!-- BGC Quality/Novelty -->
            <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin-top: 12px;">
                <div style="background: rgba(90, 138, 138, 0.15); border: 1px solid rgba(90, 138, 138, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #5a8a8a;">{kcb_stats.get('hit_percentage', 0)}%</div>
                    <div style="color: #555; font-size: 0.9em;">Match Known Clusters</div>
                    <div style="font-size: 0.8em; margin-top: 3px; color: #777;">{kcb_stats.get('regions_with_hits', 0)} of {kcb_stats.get('total_regions', 0)} BGCs</div>
                </div>
                <div style="background: rgba(106, 122, 133, 0.15); border: 1px solid rgba(106, 122, 133, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #6a7a85;">{kcb_stats.get('novel_bgc_count', 0)}</div>
                    <div style="color: #555; font-size: 0.9em;">Potentially Novel BGCs</div>
                    <div style="font-size: 0.8em; margin-top: 3px; color: #777;">No match in MIBiG</div>
                </div>
                <div style="background: rgba(122, 122, 128, 0.15); border: 1px solid rgba(122, 122, 128, 0.3); padding: 12px 15px; border-radius: 8px;">
                    <div style="font-size: 1.6em; font-weight: bold; color: #7a7a80;">{kcb_stats.get('contig_edge_count', 0)}</div>
                    <div style="color: #555; font-size: 0.9em;">Contig Edge BGCs</div>
                    <div style="font-size: 0.8em; margin-top: 3px; color: #777;">Potentially incomplete</div>
                </div>
            </div>

            {clustering_summary_html}

            <!-- BGC Type Distribution - Full Width -->
            <div class="plot" style="max-width: none;">
                <h3 style="margin-top: 0;">BGC Type Distribution</h3>
                {'<img src="bgc_donut_chart.png" alt="BGC Type Distribution" style="max-width: 100%; height: auto;">' if donut_chart_generated else '<p style="color: #999;">No BGC data available</p>'}
            </div>

            <!-- BGC Identification Status - Centered -->
            <div class="plot" style="max-width: 600px; margin: 20px auto;">
                <h3 style="margin-top: 0;">BGC Identification Status</h3>
                <img src="kcb_identification_chart.png" alt="KCB Identification Chart" onerror="this.parentElement.innerHTML='<p style=\\'color: #999;\\'>Run with --antismash_cb_knownclusters true to see identification data</p>'">
            </div>
            {kcb_mapping_section}
            {rarefaction_section}
            {phylo_tree_download_section}
        </div>

        <!-- Tab 2: Taxonomy -->
        <div class="tab-content" id="content2">
            <h2>Taxonomic Distribution of BGCs</h2>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Expandable taxonomy tree showing BGC statistics at each taxonomic level. Click on nodes to expand/collapse.
                Species nodes expand to show individual genomes with their BGC counts.</em>
            </p>
            {tree_html if tree_html else '<div class="info-box warning"><p>Taxonomy tree data not available.</p></div>'}
        </div>

        <!-- Tab 3: BGC Distribution -->
        <div class="tab-content" id="content3">
            <h2>BGC Distribution</h2>
            {tree_section}
        </div>

        <!-- Tab 4: Genomes -->
        <div class="tab-content" id="content4">
            <h2>All Genomes</h2>
            <p style="color: #666; margin-bottom: 15px;">
                <em>Searchable table of all analyzed genomes. Click genome names for detailed metadata pages.</em>
            </p>
            <div class="search-box">
                <input type="text" id="genomeSearch" placeholder="Search genomes..." onkeyup="filterGenomes()">
            </div>
            <div class="table-container">
                <table id="genomeTable">
                    <thead>
                        <tr>
                            <th>Genome Name</th>
                            <th>Assembly ID</th>
                            <th>Organism</th>
                            <th>Taxonomy</th>
                            <th>Total BGCs</th>
                            <th>Top BGC Types</th>
                        </tr>
                    </thead>
                    <tbody id="genomeTableBody">
                    {genome_table_rows}
                    </tbody>
                </table>
            </div>
        </div>

        <!-- Tab 5: Clustering -->
        <div class="tab-content" id="content5">
            <h2>BGC Clustering Analysis</h2>

            {bigscape_section_html}

            {bigslice_section_html}

            {f'<div class="info-box" style="background-color: #f8f9fa; border-left: 4px solid #6c757d;"><p style="color: #666;">No clustering analysis was performed. To enable clustering, run the pipeline with <code>--clustering bigscape</code> or <code>--clustering bigslice</code>.</p></div>' if not bigscape_section_html and not bigslice_section_html else ''}
        </div>

        <!-- Tab 6: Novel BGCs -->
        <div class="tab-content" id="content6">
            {novel_bgcs_tab_content}
        </div>

        <!-- Tab 7: KCB Hits -->
        <div class="tab-content" id="content7">
            {kcb_hits_tab_content}
        </div>

        <!-- Tab 8: Resources -->
        <div class="tab-content" id="content8">
            <h2>Pipeline Resource Usage</h2>
            <p style="color: #666; margin-bottom: 20px;">
                <em>Resource consumption metrics from Nextflow trace data, showing CPU, memory, and runtime for each pipeline process.</em>
            </p>
            {resource_usage_html if resource_usage_html else '<div class="info-box warning"><p>No resource usage data available. Trace data will appear here after running the pipeline.</p></div>'}

            {versions_html}
        </div>
    </div>

    <script>
        function filterGenomes() {{
            const input = document.getElementById('genomeSearch');
            const filter = input.value.toLowerCase();
            const tbody = document.getElementById('genomeTableBody');
            const rows = tbody.getElementsByTagName('tr');

            for (let i = 0; i < rows.length; i++) {{
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {{
                    if (cells[j].textContent.toLowerCase().indexOf(filter) > -1) {{
                        found = true;
                        break;
                    }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}

        function filterNovelBGCs() {{
            const input = document.getElementById('novelSearch');
            const filter = input.value.toLowerCase();
            const tbody = document.getElementById('novelTableBody');
            const rows = tbody.getElementsByTagName('tr');

            for (let i = 0; i < rows.length; i++) {{
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {{
                    if (cells[j].textContent.toLowerCase().indexOf(filter) > -1) {{
                        found = true;
                        break;
                    }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}

        function filterKCBHits() {{
            const input = document.getElementById('kcbSearch');
            const filter = input.value.toLowerCase();
            const tbody = document.getElementById('kcbTableBody');
            const rows = tbody.getElementsByTagName('tr');

            for (let i = 0; i < rows.length; i++) {{
                const cells = rows[i].getElementsByTagName('td');
                let found = false;
                for (let j = 0; j < cells.length; j++) {{
                    if (cells[j].textContent.toLowerCase().indexOf(filter) > -1) {{
                        found = true;
                        break;
                    }}
                }}
                rows[i].style.display = found ? '' : 'none';
            }}
        }}

        function toggleNode(nodeId) {{
            const element = document.getElementById(nodeId);
            const header = element.previousElementSibling;
            const icon = header.querySelector('.toggle-icon');
            if (element.style.display === 'none') {{
                element.style.display = 'block';
                icon.innerHTML = '&#9660;';
            }} else {{
                element.style.display = 'none';
                icon.innerHTML = '&#9654;';
            }}
        }}

        // Initialize all tree nodes as collapsed except root
        document.addEventListener('DOMContentLoaded', function() {{
            const nodeChildren = document.querySelectorAll('.node-children');
            nodeChildren.forEach(function(node, index) {{
                if (index > 0) {{
                    node.style.display = 'none';
                    const header = node.previousElementSibling;
                    if (header) {{
                        const icon = header.querySelector('.toggle-icon');
                        if (icon) icon.innerHTML = '&#9654;';
                    }}
                }}
            }});
        }});
    </script>
</body>
</html>
'''

    with open(f'{outdir}/bgc_report.html', 'w') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Visualize BGC analysis results')
    parser.add_argument('--counts', type=Path, help='Path to region_counts.tsv')
    parser.add_argument('--tabulation', type=Path, help='Path to region_tabulation.tsv')
    parser.add_argument('--assembly_info', type=Path, help='Path to assembly_info_table.txt')
    parser.add_argument('--name_map', type=Path, help='Path to name_map.json')
    parser.add_argument('--taxonomy_map', type=Path, help='Path to taxonomy_map.json')
    parser.add_argument('--taxonomy_tree', type=Path, help='Path to taxonomy_tree.json')
    parser.add_argument('--phylo_tree', type=Path, help='Path to Newick phylogenetic tree file from GTDB-Tk')
    parser.add_argument('--gtdbtk_summary', type=Path, help='Path to GTDB-Tk summary TSV file')
    parser.add_argument('--bigslice_stats', type=Path, help='Path to bigslice_statistics.json')
    parser.add_argument('--bigscape_stats', type=Path, help='Path to bigscape_statistics.json')
    parser.add_argument('--bigscape_db', type=Path, help='Path to BiG-SCAPE SQLite database for rarefaction curve')
    parser.add_argument('--gcf_data', type=Path, help='Path to gcf_representatives.json')
    parser.add_argument('--trace', type=Path, help='Path to Nextflow pipeline_trace.tsv file')
    parser.add_argument('--outdir', type=Path, required=True, help='Output directory for plots')
    parser.add_argument('--taxon', type=str, default='Unknown', help='Taxon name for report')
    parser.add_argument('--mibig_included', action='store_true', help='Whether MIBiG references were included in BiG-SCAPE analysis')
    parser.add_argument('--versions', type=Path, help='Path to software_versions.json')
    parser.add_argument('--skip_tree', action='store_true', help='Skip phylogenetic tree visualization (useful for very large datasets)')
    parser.add_argument('--outgroup', type=str, help='Outgroup taxon pattern for tree pruning (e.g., "g__Escherichia")')

    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load taxonomy data if available
    taxonomy_map_data = None
    if args.taxonomy_map and args.taxonomy_map.exists():
        print(f"Loading taxonomy map...")
        with open(args.taxonomy_map, 'r') as f:
            taxonomy_map_data = json.load(f)

    taxonomy_tree_data = None
    if args.taxonomy_tree and args.taxonomy_tree.exists():
        print(f"Loading taxonomy tree...")
        with open(args.taxonomy_tree, 'r') as f:
            taxonomy_tree_data = json.load(f)

    table_header = ''
    table_rows = ''
    stats = {}
    tree_html = ''
    genome_table_html = ''

    donut_chart_generated = False

    if args.counts:
        print(f"Generating count visualizations...")

        # Calculate summary statistics (including tabulation stats if available)
        print(f"Calculating summary statistics...")
        tabulation_path = str(args.tabulation) if args.tabulation and args.tabulation.exists() else None
        stats = calculate_summary_statistics(args.counts, tabulation_path)

        # Generate donut chart for overall BGC distribution
        print(f"Generating BGC donut chart...")
        donut_chart_generated = plot_bgc_donut_chart(args.counts, args.outdir)

        # Create genome metadata pages and genome table HTML
        if args.assembly_info and args.name_map:
            print(f"Creating genome metadata pages...")
            num_pages = create_genome_metadata_pages(args.counts, args.assembly_info, args.name_map, args.outdir, args.taxon, taxonomy_map_data, tabulation_path)
            print(f"Created {num_pages} genome metadata pages")

            # Generate searchable genome table for Genomes tab
            print(f"Generating genome table HTML...")
            genome_table_html = generate_genome_table_html(args.counts, args.assembly_info, args.name_map, taxonomy_map_data)

        # Create interactive table
        print(f"Creating interactive BGC distribution table...")
        table_header, table_rows = create_bgc_distribution_table(args.counts, args.outdir)

    # Generate taxonomy tree visualization
    if taxonomy_tree_data:
        print(f"Generating taxonomy tree visualization...")
        tree_html = generate_taxonomy_tree_html(taxonomy_tree_data)

    # Generate phylogenetic tree data for JavaScript visualization
    phylo_tree_generated = False
    phylo_tree_data = None
    gtdbtk_summary_path = str(args.gtdbtk_summary) if args.gtdbtk_summary and args.gtdbtk_summary.exists() else None
    if args.skip_tree:
        print("Skipping phylogenetic tree visualization (--skip_tree enabled)")
        phylo_tree_data = {'skipped': True, 'reason': 'user_disabled'}
    elif args.phylo_tree and args.phylo_tree.exists():
        print(f"Preparing phylogenetic tree for visualization...")
        gtdbtk_summary_path = str(args.gtdbtk_summary) if args.gtdbtk_summary and args.gtdbtk_summary.exists() else None
        counts_path = str(args.counts) if args.counts and args.counts.exists() else None

        # Prepare pruned tree data for JavaScript visualization
        phylo_tree_data = prepare_phylo_tree_for_js(
            str(args.phylo_tree),
            gtdbtk_summary_path,
            counts_path,
            str(args.outdir),
            outgroup=args.outgroup
        )

        if phylo_tree_data:
            phylo_tree_generated = True
            print(f"Prepared interactive tree with {phylo_tree_data.get('leaf_count', 0)} genomes")

    # Generate BiG-SLiCE statistics visualization
    bigslice_stats_html = ''
    if args.bigslice_stats and args.bigslice_stats.exists():
        print(f"Generating BiG-SLiCE statistics visualization...")
        bigslice_stats_html = generate_bigslice_stats_html(str(args.bigslice_stats))

    # Generate BiG-SLiCE section HTML
    bigslice_section_html = ''
    if bigslice_stats_html:
        # Sanitize taxon name to match Nextflow sanitizeTaxon function
        taxon_sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', args.taxon)
        taxon_sanitized = re.sub(r'_+', '_', taxon_sanitized).strip('_')
        bigslice_section_html = f'''<div class="clustering-section">
                <h3>BiG-SLiCE Clustering</h3>
                <div class="info-box">
                    <p>BiG-SLiCE provides scalable k-mer-based clustering with an interactive web interface.</p>
                    <p style="margin-top: 15px;"><strong>To view results:</strong></p>
                    <ol style="margin: 10px 0 10px 20px; line-height: 1.8;">
                        <li>Navigate to: <code>cd results/bigslice_results/{taxon_sanitized}/bigslice_output/</code></li>
                        <li>Start server: <code>bash start_server.sh</code></li>
                        <li>Open: <a href="http://localhost:5000" target="_blank">http://localhost:5000</a></li>
                    </ol>
                </div>
                {bigslice_stats_html}
            </div>'''

    # Generate BiG-SCAPE statistics visualization
    bigscape_stats_html = ''
    rarefaction_stats = None
    if args.bigscape_stats and args.bigscape_stats.exists():
        print(f"Generating BiG-SCAPE statistics visualization...")
        bigscape_stats_html = generate_bigscape_stats_html(str(args.bigscape_stats), args.mibig_included)

        # Generate rarefaction curve from BiG-SCAPE database
        bigscape_db = None
        if args.bigscape_db and args.bigscape_db.exists():
            bigscape_db = args.bigscape_db
        else:
            # Fallback: try to find .db file in the same directory as stats
            bigscape_dir = args.bigscape_stats.parent
            db_files = list(bigscape_dir.glob('*.db'))
            if db_files:
                bigscape_db = db_files[0]

        if bigscape_db:
            print(f"Generating rarefaction curve from {bigscape_db.name}...")
            rarefaction_stats = generate_rarefaction_curve(str(bigscape_db), args.outdir, args.taxon)
            if rarefaction_stats:
                print(f"Rarefaction curve generated: {rarefaction_stats['total_gcfs']} GCFs, {rarefaction_stats['saturation']:.0f}% saturation")
        else:
            print("Warning: BiG-SCAPE database not found, skipping rarefaction curve generation")

    # Generate GCF visualization HTML and load GCF data for overview
    gcf_visualization_html = ''
    gcf_data_dict = None
    if args.gcf_data and args.gcf_data.exists():
        print(f"Generating GCF visualization...")
        gcf_visualization_html = generate_gcf_visualization_html(str(args.gcf_data), args.taxon)
        # Also load as dict for overview sections
        try:
            with open(args.gcf_data, 'r') as f:
                gcf_data_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load GCF data for overview: {e}")

    # Load taxonomy map for overview sections
    taxonomy_map_dict = None
    if args.taxonomy_map and args.taxonomy_map.exists():
        try:
            with open(args.taxonomy_map, 'r') as f:
                taxonomy_map_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load taxonomy map: {e}")

    # Generate resource usage HTML from trace file
    resource_usage_html = ''
    if args.trace and args.trace.exists():
        print(f"Processing pipeline trace data...")
        tasks = parse_trace_file(str(args.trace))
        if tasks:
            processes = aggregate_trace_by_process(tasks)
            resource_usage_html = generate_resource_usage_html(tasks, processes)
            print(f"Processed {len(tasks)} tasks from trace file")

    if args.tabulation and args.tabulation.exists():
        print(f"Generating KCB identification chart...")
        plot_kcb_identification_chart(stats.get('kcb_stats', {}), args.outdir)

    # Load software versions if available
    versions_data = None
    if args.versions and args.versions.exists():
        print(f"Loading software versions...")
        try:
            with open(args.versions, 'r') as f:
                versions_data = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load versions file: {e}")

    # Load name_map for BGC Distribution tab
    name_map_dict = None
    if args.name_map and args.name_map.exists():
        try:
            with open(args.name_map, 'r') as f:
                name_map_dict = json.load(f)
        except Exception as e:
            print(f"Warning: Could not load name_map file: {e}")

    if args.counts or args.tabulation:
        print(f"Generating HTML report...")
        generate_html_report(args.outdir, args.taxon, table_header, table_rows, stats, tree_html,
                            bigslice_stats_html, bigscape_stats_html, gcf_visualization_html,
                            donut_chart_generated, phylo_tree_generated,
                            genome_table_html, resource_usage_html, phylo_tree_data,
                            gcf_data=gcf_data_dict, taxonomy_map=taxonomy_map_dict,
                            bigslice_section_html=bigslice_section_html,
                            versions_data=versions_data,
                            rarefaction_stats=rarefaction_stats,
                            gtdbtk_summary_path=gtdbtk_summary_path,
                            name_map=name_map_dict)
        print(f"Visualizations complete! Open {args.outdir}/bgc_report.html in a browser.")

if __name__ == '__main__':
    main()
