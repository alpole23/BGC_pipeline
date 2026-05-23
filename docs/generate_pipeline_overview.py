#!/usr/bin/env python3
"""Generate the pipeline overview figure for ClusterQuest."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

fig, ax = plt.subplots(1, 1, figsize=(16, 14))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Color scheme
colors = {
    'genome': '#6baed6',      # Blue - Genome Processing
    'bgc': '#e6854a',         # Orange - BGC Detection
    'clustering': '#74c476',  # Green - GCF Clustering
    'phylogeny': '#9e7fc0',   # Purple - Phylogenetics
    'viz': '#7a9eb8',         # Slate blue - Visualization
    'output': '#f0c75e',      # Yellow-gold - Output files
    'database': '#bdbdbd',    # Gray - Databases
}

# Title
ax.text(50, 97, 'ClusterQuest Pipeline', fontsize=24, fontweight='bold',
        ha='center', va='top')
ax.text(50, 93, 'Large-scale Organization Of Secondary Metabolites', fontsize=14,
        ha='center', va='top', style='italic')

# Helper function for process boxes
def add_process_box(ax, x, y, width, height, color, title, subtitle):
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.5",
                         facecolor=color, edgecolor='#333333', linewidth=2)
    ax.add_patch(box)
    ax.text(x, y + 2, title, fontsize=14, fontweight='bold',
            ha='center', va='center')
    ax.text(x, y - 3, subtitle, fontsize=11, ha='center', va='center',
            style='italic', color='#333333')

# Helper function for database cylinders
def add_database(ax, x, y, width, height, label):
    # Draw cylinder body
    from matplotlib.patches import Ellipse, Rectangle
    body = Rectangle((x - width/2, y - height/2), width, height,
                     facecolor=colors['database'], edgecolor='#333333', linewidth=1.5)
    ax.add_patch(body)
    # Top ellipse
    top = Ellipse((x, y + height/2), width, height*0.3,
                  facecolor=colors['database'], edgecolor='#333333', linewidth=1.5)
    ax.add_patch(top)
    # Bottom ellipse (partial)
    bottom = Ellipse((x, y - height/2), width, height*0.3,
                     facecolor=colors['database'], edgecolor='#333333', linewidth=1.5)
    ax.add_patch(bottom)
    ax.text(x, y, label, fontsize=9, ha='center', va='center', fontweight='bold')

# Helper function for output boxes
def add_output_box(ax, x, y, width, height, label):
    box = FancyBboxPatch((x - width/2, y - height/2), width, height,
                         boxstyle="round,pad=0.02,rounding_size=0.3",
                         facecolor=colors['output'], edgecolor='#333333', linewidth=1.5)
    ax.add_patch(box)
    ax.text(x, y, label, fontsize=11, ha='center', va='center', color='#f5f5f0', fontweight='bold')

# Helper function for arrows
def add_arrow(ax, start, end, color='#2d8a2d', style='-', linewidth=2):
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='->', color=color, lw=linewidth,
                               linestyle=style, shrinkA=0, shrinkB=0))

# Process positions (x, y)
positions = {
    'download': (50, 82),
    'antismash': (50, 66),
    'clustering': (50, 50),
    'phylogeny': (50, 34),
    'visualize': (50, 18),
}

# Add process boxes
add_process_box(ax, *positions['download'], 30, 10, colors['genome'],
                'DOWNLOAD_GENOMES', 'NCBI Datasets + Taxonomy')
add_process_box(ax, *positions['antismash'], 30, 10, colors['bgc'],
                'ANTISMASH_ANALYSIS', 'Phosphonate BGC Detection')
add_process_box(ax, *positions['clustering'], 30, 10, colors['clustering'],
                'CLUSTERING', 'BiG-SCAPE GCF Clustering')
add_process_box(ax, *positions['phylogeny'], 30, 10, colors['phylogeny'],
                'PHYLOGENY', 'GTDB-Tk Classification')
add_process_box(ax, *positions['visualize'], 30, 10, colors['viz'],
                'VISUALIZE_RESULTS', 'Interactive HTML Report')

# Add databases (left side)
add_database(ax, 12, 82, 12, 8, 'NCBI')
add_database(ax, 12, 66, 12, 8, 'antiSMASH DB')
add_database(ax, 12, 50, 12, 8, 'Pfam / MIBiG')
add_database(ax, 12, 34, 12, 8, 'GTDB-Tk DB')

# Add output boxes (right side)
add_output_box(ax, 88, 82, 18, 6, 'Renamed Genomes')
add_output_box(ax, 88, 66, 18, 6, 'BGC Results (.json)')
add_output_box(ax, 88, 50, 18, 6, 'GCF Clusters (.db)')
add_output_box(ax, 88, 34, 18, 6, 'Phylo Tree (.nwk)')
add_output_box(ax, 88, 18, 18, 6, 'Report (.html)')

# Main flow arrows (vertical, green)
add_arrow(ax, (50, 77), (50, 71))
add_arrow(ax, (50, 61), (50, 55))
add_arrow(ax, (50, 45), (50, 39))
add_arrow(ax, (50, 29), (50, 23))

# Database to process arrows (horizontal, green)
add_arrow(ax, (18, 82), (35, 82))
add_arrow(ax, (18, 66), (35, 66))
add_arrow(ax, (18, 50), (35, 50))
add_arrow(ax, (18, 34), (35, 34))

# Process to output arrows (horizontal, green)
add_arrow(ax, (65, 82), (79, 82))
add_arrow(ax, (65, 66), (79, 66))
add_arrow(ax, (65, 50), (79, 50))
add_arrow(ax, (65, 34), (79, 34))
add_arrow(ax, (65, 18), (79, 18))

# Cross-taxon reuse arrows (red dashed) - from output boxes to process boxes
# reuse_antismash_from: BGC Results -> ANTISMASH_ANALYSIS
ax.annotate('', xy=(65, 63), xytext=(79, 63),
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2.5,
                           linestyle='--', connectionstyle='arc3,rad=-0.3'))

# reuse_gtdbtk_from: Phylo Tree -> PHYLOGENY
ax.annotate('', xy=(65, 31), xytext=(79, 31),
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2.5,
                           linestyle='--', connectionstyle='arc3,rad=-0.3'))

# Legend - using Line2D for proper line style representation
legend_y = 5
legend_elements = [
    mpatches.Patch(facecolor=colors['genome'], edgecolor='#333333', label='Genome\nProcessing'),
    mpatches.Patch(facecolor=colors['bgc'], edgecolor='#333333', label='BGC\nDetection'),
    mpatches.Patch(facecolor=colors['clustering'], edgecolor='#333333', label='GCF\nClustering'),
    mpatches.Patch(facecolor=colors['phylogeny'], edgecolor='#333333', label='Phylogenetics'),
    mpatches.Patch(facecolor=colors['viz'], edgecolor='#333333', label='Visualization'),
    Line2D([0], [0], color='#c0392b', linestyle='--', linewidth=2,
           label='Cross-taxon\nresult reuse'),
]

legend = ax.legend(handles=legend_elements, loc='lower center', ncol=6,
                   fontsize=11, frameon=False, bbox_to_anchor=(0.5, 0.01),
                   handlelength=2.5, handleheight=1.5)

plt.savefig('docs/pipeline_overview.png', dpi=130, facecolor='white',
            bbox_inches='tight', pad_inches=0.3)
print("Pipeline overview figure saved to docs/pipeline_overview.png")
