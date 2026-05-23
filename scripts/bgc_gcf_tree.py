#!/usr/bin/env python3
"""
GCF biosynthetic phylogeny using center-to-center BiG-SCAPE distances.

Extracts the representative (center) BGC for each GCF from the BiG-SCAPE
SQLite database, looks up their pairwise distances, and builds a
Neighbor-Joining tree. The center BGC is the medoid — the member with the
minimum sum of distances to all other GCF members — as defined by BiG-SCAPE.

Usage:
    python scripts/bgc_gcf_tree.py \\
        --db results/bigscape_results/Pantoea/Pantoea.db \\
        --outdir results/bgc_trees/Pantoea \\
        [--cutoff 0.3] \\
        [--coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt]

Outputs:
    gcf_biosynthetic_tree.nwk        — NJ tree in Newick format
    gcf_center_distances.tsv         — full pairwise distance matrix
    gcf_biosynthetic_tree.png/.svg   — figure (when --coupling_annotation provided)
"""

import argparse
import csv
import os
import re
import sqlite3
from collections import Counter, defaultdict
from itertools import combinations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import COUPLING_COLORS, COUPLING_ORDER


def get_gcf_info(conn, cutoff):
    """Return {gcf_id: {center_id, n_members}} for region-level GCFs."""
    cur = conn.cursor()
    cur.execute("""
        SELECT f.id, f.center_id, COUNT(rf.record_id) AS n
        FROM family f
        JOIN bgc_record_family rf ON rf.family_id = f.id
        JOIN bgc_record br        ON br.id = rf.record_id
        WHERE f.cutoff = ? AND br.record_type = 'region'
        GROUP BY f.id
        ORDER BY n DESC
    """, (cutoff,))
    return {row[0]: {'center_id': row[1], 'n_members': row[2]}
            for row in cur.fetchall()}


def get_center_distances(conn, gcf_info):
    """
    Look up the BiG-SCAPE pairwise distance between every pair of GCF centers.
    Returns {(gcf_a, gcf_b): distance} for gcf_a < gcf_b.
    """
    cur = conn.cursor()
    gcf_ids = sorted(gcf_info)
    distances = {}
    missing = []

    for gcf_a, gcf_b in combinations(gcf_ids, 2):
        ca = gcf_info[gcf_a]['center_id']
        cb = gcf_info[gcf_b]['center_id']
        cur.execute("""
            SELECT distance FROM distance
            WHERE (record_a_id = ? AND record_b_id = ?)
               OR (record_a_id = ? AND record_b_id = ?)
            LIMIT 1
        """, (ca, cb, cb, ca))
        row = cur.fetchone()
        if row:
            distances[(gcf_a, gcf_b)] = row[0]
        else:
            missing.append((gcf_a, gcf_b))
            distances[(gcf_a, gcf_b)] = 1.0   # maximum distance as fallback

    if missing:
        print(f'  Warning: {len(missing)} center pairs had no stored distance '
              f'(set to 1.0): {missing}')
    return distances


def build_nj_tree(gcf_ids, distances):
    """Build a Neighbor-Joining tree from the center-to-center distance matrix."""
    labels = [f'GCF-{gid}' for gid in gcf_ids]

    # Biopython DistanceMatrix expects a lower-triangular list of lists
    dm_rows = []
    for i, gcf_a in enumerate(gcf_ids):
        row = []
        for j in range(i + 1):
            gcf_b = gcf_ids[j]
            if gcf_a == gcf_b:
                row.append(0.0)
            else:
                key = (min(gcf_a, gcf_b), max(gcf_a, gcf_b))
                row.append(distances[key])
        dm_rows.append(row)

    dm = DistanceMatrix(labels, dm_rows)
    tree = DistanceTreeConstructor().nj(dm)
    return tree, dm


def save_distance_matrix(gcf_ids, distances, path):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([''] + [f'GCF-{gid}' for gid in gcf_ids])
        for gcf_a in gcf_ids:
            row = [f'GCF-{gcf_a}']
            for gcf_b in gcf_ids:
                if gcf_a == gcf_b:
                    row.append('0.000000')
                else:
                    key = (min(gcf_a, gcf_b), max(gcf_a, gcf_b))
                    row.append(f'{distances[key]:.6f}')
            writer.writerow(row)


def load_coupling_classes(path):
    """Parse iTOL coupling annotation → {gbk_basename: class_name}."""
    classes = {}
    legacy = {'Fe-ADH': 'VlpB-like', 'TPP+NTP': 'Ppd-CDP',
               'PalB': 'PalB-like', 'FrbC': 'FrbC-like'}
    in_data = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == 'DATA':
                in_data = True; continue
            if in_data and line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    label = parts[0]
                    cls   = legacy.get(parts[2], parts[2])
                    if not any(label.endswith(f'_{i}') for i in range(10)):
                        classes[label] = cls
    return classes


def get_gcf_coupling_class(conn, gcf_info, cutoff, coupling_classes):
    """Return dominant coupling class for each GCF using member BGC labels."""
    cur = conn.cursor()
    gcf_dominant = {}
    for gcf_id, info in gcf_info.items():
        cur.execute("""
            SELECT g.path
            FROM bgc_record_family rf
            JOIN bgc_record br ON br.id = rf.record_id
            JOIN gbk g ON g.id = br.gbk_id
            WHERE rf.family_id = ? AND br.record_type = 'region'
        """, (gcf_id,))
        counts = Counter()
        for (path,) in cur.fetchall():
            gbk_base = os.path.splitext(os.path.basename(path))[0]
            counts[coupling_classes.get(gbk_base, 'Unknown')] += 1
        gcf_dominant[gcf_id] = counts.most_common(1)[0][0] if counts else 'Unknown'
    return gcf_dominant


def assign_layout(clade, counter, depth=0):
    clade._depth = depth
    if clade.is_terminal():
        clade._x = counter[0]
        counter[0] += 1
        return
    for child in clade.clades:
        assign_layout(child, counter, depth + 1)
    clade._x = sum(c._x for c in clade.clades) / len(clade.clades)


def max_depth(clade):
    if clade.is_terminal():
        return clade._depth
    return max(max_depth(c) for c in clade.clades)


def draw_cladogram(ax, clade, color='#333333', lw=1.2):
    if clade.is_terminal():
        return
    x_node   = clade._depth
    child_ys = [c._x for c in clade.clades]
    ax.plot([x_node, x_node], [min(child_ys), max(child_ys)],
            color=color, lw=lw, solid_capstyle='round')
    for child in clade.clades:
        ax.plot([x_node, child._depth], [child._x, child._x],
                color=color, lw=lw, solid_capstyle='round')
        draw_cladogram(ax, child, color, lw)


def plot_gcf_tree(tree, gcf_info, gcf_dominant, outdir):
    """Draw the GCF biosynthetic NJ tree with coupling-enzyme coloring."""
    assign_layout(tree.root, [0])
    md = max_depth(tree.root)

    leaf_order = list(tree.get_terminals())
    n = len(leaf_order)

    fig_w = 6.0
    fig_h = max(3.0, n * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    draw_cladogram(ax, tree.root)

    # Draw leaf nodes as colored circles + labels
    max_members = max(info['n_members'] for info in gcf_info.values())
    for clade in leaf_order:
        name = clade.name.strip("'") if clade.name else ''
        try:
            gcf_id = int(name.replace('GCF-', ''))
        except ValueError:
            gcf_id = None

        y     = clade._x
        x     = clade._depth
        cls   = gcf_dominant.get(gcf_id, 'Unknown')
        color = COUPLING_COLORS.get(cls, '#aaaaaa')
        n_mem = gcf_info[gcf_id]['n_members'] if gcf_id else 1

        # Circle size proportional to member count (area ∝ n)
        size = 40 + 300 * (n_mem / max_members)
        ax.scatter(x, y, s=size, color=color, zorder=5,
                   edgecolors='white', linewidths=0.8)

        label = f'GCF-{gcf_id}  (n={n_mem})' if gcf_id else name
        ax.text(x + 0.15, y, label, va='center', ha='left',
                fontsize=9, color='#222222')

    ax.set_xlim(-0.3, md + 2.5)
    ax.set_ylim(-0.5, n - 0.5)
    ax.axis('off')

    # Combined legend: coupling enzyme classes + GCF size scale
    present_classes = set(gcf_dominant.values())
    patches = [mpatches.Patch(color=COUPLING_COLORS[c], label=c)
               for c in COUPLING_ORDER if c in present_classes]
    size_handles = [
        plt.scatter([], [], s=40 + 300*(sz/max_members), color='#888',
                    label=f'n={sz}', edgecolors='white', linewidths=0.8)
        for sz in [1, max(max_members // 2, 1), max_members]
    ]
    ax.legend(handles=patches + size_handles,
              title='Coupling enzyme / GCF size',
              title_fontsize=8, fontsize=8, framealpha=0.9,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              bbox_transform=ax.transAxes, borderaxespad=0)

    fig.suptitle('GCF Biosynthetic Phylogeny\n'
                 'NJ tree · BiG-SCAPE center-to-center distances',
                 fontsize=10, y=0.98, va='top')

    out_png = os.path.join(outdir, 'gcf_biosynthetic_tree.png')
    out_svg = os.path.join(outdir, 'gcf_biosynthetic_tree.svg')
    fig.savefig(out_png, dpi=180, bbox_inches='tight')
    fig.savefig(out_svg,           bbox_inches='tight')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='NJ tree of GCFs from BiG-SCAPE center-to-center distances',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--db',      required=True, help='BiG-SCAPE SQLite database')
    parser.add_argument('--outdir',  required=True, help='Output directory')
    parser.add_argument('--cutoff',  type=float, default=0.3,
                        help='GCF cutoff used in BiG-SCAPE run (default: 0.3)')
    parser.add_argument('--coupling_annotation', default=None,
                        help='iTOL coupling annotation file for figure coloring')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    conn = sqlite3.connect(args.db)

    print('Loading GCF info...')
    gcf_info = get_gcf_info(conn, args.cutoff)
    gcf_ids  = sorted(gcf_info)
    print(f'  {len(gcf_ids)} GCFs at cutoff {args.cutoff}:')
    for gid in gcf_ids:
        info = gcf_info[gid]
        print(f'    GCF-{gid}: n={info["n_members"]:4d}  center_bgc_id={info["center_id"]}')

    print('\nLooking up center-to-center distances...')
    distances = get_center_distances(conn, gcf_info)

    print('\nDistance matrix (center-to-center):')
    header = '       ' + '  '.join(f'GCF-{g:<2}' for g in gcf_ids)
    print(header)
    for gcf_a in gcf_ids:
        row_vals = []
        for gcf_b in gcf_ids:
            if gcf_a == gcf_b:
                row_vals.append('  0.000')
            else:
                key = (min(gcf_a, gcf_b), max(gcf_a, gcf_b))
                row_vals.append(f'  {distances[key]:.3f}')
        print(f'GCF-{gcf_a:<2}' + ''.join(row_vals))

    out_tsv = os.path.join(args.outdir, 'gcf_center_distances.tsv')
    save_distance_matrix(gcf_ids, distances, out_tsv)
    print(f'\nSaved distance matrix: {out_tsv}')

    print('\nBuilding Neighbor-Joining tree...')
    tree, dm = build_nj_tree(gcf_ids, distances)

    out_nwk = os.path.join(args.outdir, 'gcf_biosynthetic_tree.nwk')
    Phylo.write(tree, out_nwk, 'newick')
    print(f'Saved Newick tree:     {out_nwk}')

    print('\nTree (ASCII):')
    Phylo.draw_ascii(tree)

    if args.coupling_annotation:
        print('\nGenerating figure...')
        coupling_classes = load_coupling_classes(args.coupling_annotation)
        gcf_dominant = get_gcf_coupling_class(conn, gcf_info, args.cutoff, coupling_classes)
        plot_gcf_tree(tree, gcf_info, gcf_dominant, args.outdir)

    conn.close()
    print('Done.')


if __name__ == '__main__':
    main()
