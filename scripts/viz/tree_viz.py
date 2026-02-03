#!/usr/bin/env python3
"""Tree visualization functions for phylogenetic and taxonomic trees."""

import re
import sys
import traceback
from collections import defaultdict
from io import StringIO
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from viz.charts import get_bgc_color


# =============================================================================
# NEWICK PARSING UTILITIES
# =============================================================================

def parse_newick(newick_str):
    """
    Parse a Newick format string into a tree structure.
    Returns a dictionary representing the tree with nodes containing:
    - name: node name (leaf names or internal node labels)
    - branch_length: distance from parent
    - children: list of child nodes
    """
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
    Returns list of (leaf_name, cumulative_distance, node) tuples.
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


# =============================================================================
# TAXONOMY TREE VISUALIZATION
# =============================================================================

def plot_circular_taxonomy_tree(taxonomy_tree_data, counts_file, outdir):
    """Create circular taxonomy tree with BGC annotation rings."""
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
        """Recursively collect leaf nodes with their taxonomic paths"""
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

def plot_circular_phylogenetic_tree(newick_file, counts_file, gtdbtk_summary, outdir):
    """
    Create circular phylogenetic tree from GTDB-Tk Newick output with BGC annotation rings.

    Args:
        newick_file: Path to Newick tree file from GTDB-Tk
        counts_file: Path to region_counts.tsv for BGC data
        gtdbtk_summary: Path to GTDB-Tk summary TSV for genome name mapping
        outdir: Output directory for the plot
    """
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


def prepare_phylo_tree_for_js(newick_file, gtdbtk_summary, counts_file, outdir):
    """
    Prepare phylogenetic tree data for JavaScript visualization.
    Prunes the GTDB-Tk tree to only include user genomes and exports data for JS rendering.
    Uses Bio.Phylo for efficient parsing of large trees.

    Args:
        newick_file: Path to Newick tree file from GTDB-Tk
        gtdbtk_summary: Path to GTDB-Tk summary TSV
        counts_file: Path to region_counts.tsv for BGC data
        outdir: Output directory

    Returns:
        dict with 'newick' (pruned tree string) and 'metadata' (genome info), or None on failure
    """
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
        from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

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

        # Prune tree to user genomes only
        print(f"Pruning tree to {len(user_terminals)} user genomes...")

        # For small numbers of genomes, build a simple tree with just their relationships
        if len(user_terminals) == 1:
            # Single genome - create simple tree
            name = user_terminals[0].name
            bl = user_terminals[0].branch_length or 0.0
            pruned_newick = f"({name}:{bl});"
        else:
            # Build tree from pairwise distances using neighbor-joining
            # Get pairwise distances using tree structure
            names = [t.name for t in user_terminals]
            n = len(names)
            print(f"Computing pairwise distances for {n} genomes...")

            # Build distance matrix from tree distances
            # For large trees, compute in batches to show progress
            matrix = []
            for i in range(n):
                if i % 50 == 0 and i > 0:
                    print(f"  Distance matrix: {i}/{n} rows computed...")
                row = []
                for j in range(i + 1):
                    if i == j:
                        row.append(0.0)
                    else:
                        # Get distance between two terminals in tree
                        try:
                            dist = tree.distance(user_terminals[i], user_terminals[j])
                        except:
                            dist = 1.0  # Default if distance calc fails
                        row.append(dist)
                matrix.append(row)

            print(f"Building neighbor-joining tree from distance matrix...")
            # Create distance matrix and build tree
            dm = DistanceMatrix(names, matrix)
            constructor = DistanceTreeConstructor()
            pruned_tree_obj = constructor.nj(dm)

            # Fix negative branch lengths (common in NJ trees)
            for clade in pruned_tree_obj.find_clades():
                if clade.branch_length is not None and clade.branch_length < 0:
                    clade.branch_length = 0.0001  # Small positive value

            # Write to Newick
            output = StringIO()
            Phylo.write(pruned_tree_obj, output, 'newick')
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
    leaf_pattern = re.compile(r'[(),]([A-Za-z_][^:(),]*):')
    estimated_leaves = len(leaf_pattern.findall(newick_str))
    print(f"Estimated tree size: ~{estimated_leaves} leaves")

    # Skip if too large for custom parser
    MAX_TREE_SIZE = 5000
    if estimated_leaves > MAX_TREE_SIZE:
        print(f"Tree too large for custom parser ({estimated_leaves} > {MAX_TREE_SIZE})")
        print("Install biopython for large tree support: pip install biopython")
        return {'skipped': True, 'reason': 'tree_too_large', 'leaf_count': estimated_leaves, 'user_genome_count': len(user_genomes)}

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


def generate_static_circular_tree(newick_str, outdir, metadata=None):
    """
    Generate a static circular/radial phylogenetic tree image using matplotlib.

    Args:
        newick_str: Newick format tree string
        outdir: Output directory for the image
        metadata: Optional dict of genome metadata

    Returns:
        Path to generated image, or None on failure
    """
    try:
        from Bio import Phylo

        # Parse tree
        tree = Phylo.read(StringIO(newick_str), 'newick')
        terminals = tree.get_terminals()
        n_leaves = len(terminals)

        if n_leaves == 0:
            print("No leaves in tree for static image")
            return None

        print(f"Generating static circular tree image with {n_leaves} leaves...")

        # Calculate figure size based on number of leaves
        # More leaves = larger figure for readability
        if n_leaves <= 20:
            figsize = (10, 10)
            fontsize = 9
        elif n_leaves <= 50:
            figsize = (12, 12)
            fontsize = 7
        elif n_leaves <= 150:
            figsize = (16, 16)
            fontsize = 5
        else:
            figsize = (20, 20)
            fontsize = 3

        # Create figure with polar projection for circular tree
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'polar'}, facecolor='white')
        ax.set_facecolor('white')

        # Get tree depth for scaling
        depths = tree.depths()
        max_depth = max(depths.values()) if depths else 1.0
        if max_depth == 0:
            max_depth = 1.0

        # Assign angular positions to leaves
        leaf_angles = {}
        for i, terminal in enumerate(terminals):
            angle = 2 * np.pi * i / n_leaves
            leaf_angles[terminal] = angle

        # Calculate internal node angles (average of descendants)
        def get_clade_angle(clade):
            if clade.is_terminal():
                return leaf_angles[clade]
            child_angles = [get_clade_angle(c) for c in clade.clades]
            # Handle angle wraparound
            if max(child_angles) - min(child_angles) > np.pi:
                child_angles = [a if a > np.pi else a + 2*np.pi for a in child_angles]
            avg = np.mean(child_angles)
            if avg > 2*np.pi:
                avg -= 2*np.pi
            return avg

        # Draw tree branches
        def draw_clade(clade, parent_radius=0, parent_angle=0):
            clade_angle = get_clade_angle(clade)
            clade_depth = depths.get(clade, 0)
            clade_radius = (clade_depth / max_depth) * 0.8  # Scale to 80% of radius

            # Draw arc from parent angle to this clade's angle at parent's radius
            if parent_radius > 0:
                # Arc at parent radius
                if abs(clade_angle - parent_angle) > 0.01:
                    arc_angles = np.linspace(min(parent_angle, clade_angle),
                                            max(parent_angle, clade_angle), 20)
                    ax.plot(arc_angles, [parent_radius]*len(arc_angles),
                           color='#2c5aa0', linewidth=0.8, alpha=0.7)

            # Draw radial line from parent radius to this clade's radius
            ax.plot([clade_angle, clade_angle], [parent_radius, clade_radius],
                   color='#2c5aa0', linewidth=0.8, alpha=0.7)

            # Draw children
            for child in clade.clades:
                draw_clade(child, clade_radius, clade_angle)

            # Draw leaf labels
            if clade.is_terminal():
                # Add small dot at leaf position
                ax.scatter([clade_angle], [clade_radius], s=15, c='#2c5aa0', zorder=5)

                # Add label
                label = clade.name.replace('_', ' ') if clade.name else ''
                if len(label) > 30:
                    label = label[:27] + '...'

                # Rotate label to be readable
                rotation = np.degrees(clade_angle)
                if 90 < rotation < 270:
                    rotation += 180
                    ha = 'right'
                else:
                    ha = 'left'

                ax.text(clade_angle, clade_radius + 0.05, label,
                       fontsize=fontsize, rotation=rotation, rotation_mode='anchor',
                       ha=ha, va='center', color='#333')

        # Draw the tree starting from root
        draw_clade(tree.root)

        # Style the plot
        ax.set_ylim(0, 1.1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)

        # Add title
        ax.set_title(f'Phylogenetic Tree ({n_leaves} genomes)',
                    fontsize=12, color='#333', pad=20)

        # Save figure
        output_path = Path(outdir) / 'circular_phylo_tree.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"Saved static circular tree to {output_path}")
        return str(output_path)

    except ImportError as e:
        print(f"Could not generate static tree image: {e}")
        return None
    except Exception as e:
        print(f"Error generating static tree image: {e}")
        traceback.print_exc()
        return None
