#!/usr/bin/env python3
"""Chart generation functions for BGC visualization."""

import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.constants import BGC_COLORS


def get_bgc_color(bgc_type):
    """Get color for a BGC type, handling hybrids and unknown types."""
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


def plot_bgc_donut_chart(counts_file, outdir):
    """Create donut chart showing overall BGC type distribution with consistent colors."""
    df = pd.read_csv(counts_file, sep='\t', comment='#')

    # Select BGC columns
    numeric_cols = df.select_dtypes(include='number').columns
    bgc_cols = [col for col in numeric_cols if col not in ['total_count']]

    # Sum total counts for each BGC type
    bgc_totals = df[bgc_cols].sum().sort_values(ascending=False)

    # Filter out zero values
    bgc_totals = bgc_totals[bgc_totals > 0]

    if len(bgc_totals) == 0:
        print("No BGC data for donut chart")
        return False

    # Get colors for each BGC type
    colors = [get_bgc_color(bgc_type) for bgc_type in bgc_totals.index]

    # Create figure - same size as KCB chart for consistency
    fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
    ax.set_facecolor('white')

    # Create donut chart
    wedges, texts, autotexts = ax.pie(
        bgc_totals.values,
        labels=None,  # We'll add a legend instead
        colors=colors,
        autopct=lambda pct: f'{pct:.1f}%' if pct > 3 else '',
        pctdistance=0.75,
        startangle=90,
        wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2)
    )

    # Style the percentage labels
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(14)

    # Add center text with total
    total_bgcs = bgc_totals.sum()
    ax.text(0, 0, f'{int(total_bgcs)}\nTotal BGCs',
            ha='center', va='center', fontsize=26, fontweight='bold', color='#333')

    # Create legend with counts
    legend_labels = [f'{bgc_type}: {int(count)} ({count/total_bgcs*100:.1f}%)'
                     for bgc_type, count in bgc_totals.items()]
    ax.legend(wedges, legend_labels, title='BGC Types', loc='center left',
              bbox_to_anchor=(1.0, 0.5), fontsize=16, title_fontsize=18)

    plt.tight_layout()
    plt.savefig(f'{outdir}/bgc_donut_chart.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    return True


def plot_kcb_identification_chart(kcb_stats, outdir):
    """Create pie chart showing proportion of identified vs unidentified BGCs."""
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
                legend_text += f'  - {sim_level.capitalize()}: {sim_breakdown[sim_level]}\n'
        ax.annotate(legend_text.strip(), xy=(0.5, -0.12), xycoords='axes fraction',
                    ha='center', va='top', fontsize=12, color='#666',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8f9fa', edgecolor='#ddd'))

    plt.tight_layout()
    plt.savefig(f'{outdir}/kcb_identification_chart.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
