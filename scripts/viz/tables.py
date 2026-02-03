#!/usr/bin/env python3
"""Table generation functions for BGC visualization."""

import json
import sys
from pathlib import Path

import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def get_genome_count(counts_file):
    """Get the number of genomes in the counts file."""
    counts_df = pd.read_csv(counts_file, sep='\t', comment='#')
    return len(counts_df)


def generate_genome_table_html(counts_file, assembly_info, name_map, taxonomy_map_data=None):
    """Generate HTML for a searchable genome table."""
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
    """Calculate summary statistics for the dataset including tabulation stats."""
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
    """Create interactive HTML table for BGC distribution with clickable genome links and color-coding."""
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
        """Generate color based on value using sequential green palette."""
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
