#!/usr/bin/env python3
"""
Standalone NJ tree of ALL phosphonate BGCs using BiG-SCAPE pairwise distances.

Leaves (one per region-level BGC) are annotated with two colored strips:
  - Coupling enzyme class (from iTOL coupling annotation file)
  - GCF membership (from BiG-SCAPE DB)

BiG-SCAPE stores all pairwise distances (not just those below the cutoff),
so the full distance matrix is available for all 305 BGCs.

Usage:
    python scripts/bgc_all_bgcs_tree.py \\
        --db      work/57/.../Pantoea.db \\
        --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \\
        --outdir  results/bgc_trees/Pantoea \\
        [--cutoff 0.3]
"""

import argparse
import os
import sqlite3
import sys
from collections import Counter
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import COUPLING_COLORS, COUPLING_ORDER

SINGLETON_COLOR = '#dddddd'  # kept for potential future use


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_coupling_classes(path):
    """Parse iTOL coupling annotation → {gbk_basename: class_name}."""
    classes = {}
    in_data = False
    legacy = {'Fe-ADH': 'VlpB-like', 'TPP+NTP': 'Ppd-CDP',
               'PalB': 'PalB-like', 'FrbC': 'FrbC-like'}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == 'DATA':
                in_data = True
                continue
            if in_data and line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    label = parts[0]
                    cls   = legacy.get(parts[2], parts[2])
                    # Only keep region-level labels (no _1 _2 _3 suffixes)
                    if not any(label.endswith(f'_{i}') for i in range(10)):
                        classes[label] = cls
    return classes


def load_bgc_data(conn, cutoff):
    """
    Return:
      bgc_ids   : sorted list of region-level bgc_record ids
      id_to_meta: {record_id: {genome, gbk_basename, gcf_id, organism}}
      distances : {(id_a, id_b): distance}  id_a < id_b
    """
    cur = conn.cursor()

    # All region-level records with GCF assignment
    cur.execute("""
        SELECT br.id, g.path, g.organism, rf.family_id
        FROM bgc_record br
        JOIN gbk g ON g.id = br.gbk_id
        LEFT JOIN bgc_record_family rf ON rf.record_id = br.id
        WHERE br.record_type = 'region'
        ORDER BY br.id
    """)
    rows = cur.fetchall()

    id_to_meta = {}
    for rec_id, path, organism, family_id in rows:
        genome      = os.path.basename(os.path.dirname(path))
        gbk_base    = os.path.splitext(os.path.basename(path))[0]
        id_to_meta[rec_id] = {
            'genome':      genome,
            'gbk_basename': gbk_base,
            'gcf_id':      family_id,   # None for singletons
            'organism':    organism or '',
        }

    bgc_ids = sorted(id_to_meta)

    # All pairwise distances between region-level BGCs
    id_set = set(bgc_ids)
    cur.execute("""
        SELECT d.record_a_id, d.record_b_id, d.distance
        FROM distance d
        JOIN bgc_record br_a ON br_a.id = d.record_a_id
        JOIN bgc_record br_b ON br_b.id = d.record_b_id
        WHERE br_a.record_type = 'region' AND br_b.record_type = 'region'
    """)
    distances = {}
    for a, b, dist in cur.fetchall():
        key = (min(a, b), max(a, b))
        distances[key] = dist

    return bgc_ids, id_to_meta, distances


# ─── Tree building ─────────────────────────────────────────────────────────────

def build_nj_tree(bgc_ids, distances):
    """Build NJ tree from full pairwise distance matrix."""
    n      = len(bgc_ids)
    labels = [str(i) for i in bgc_ids]   # use record IDs as leaf names

    dm_rows = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                a, b = bgc_ids[j], bgc_ids[i]
                key  = (min(a, b), max(a, b))
                row.append(distances.get(key, 1.0))
        dm_rows.append(row)

    print(f'  Building NJ tree for {n} BGCs...')
    dm   = DistanceMatrix(labels, dm_rows)
    tree = DistanceTreeConstructor().nj(dm)
    return tree


# ─── Layout helpers ───────────────────────────────────────────────────────────

def assign_layout(clade, counter, depth=0):
    """Assign ._x (leaf y-position) and ._depth to every node."""
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


def draw_cladogram(ax, clade, color='#333333', lw=0.6):
    """Root on left, leaves on right."""
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


# ─── Figure ───────────────────────────────────────────────────────────────────

def plot_tree(tree, bgc_ids, id_to_meta, coupling_classes, outdir):
    # Assign layout
    assign_layout(tree.root, [0])
    md = max_depth(tree.root)

    # Leaf order from tree traversal
    leaf_order = [int(clade.name) for clade in tree.get_terminals()]
    n_leaves   = len(leaf_order)

    # Unique GCF IDs (sorted by size desc, singletons last)
    gcf_counts  = Counter(id_to_meta[i]['gcf_id'] for i in bgc_ids
                          if id_to_meta[i]['gcf_id'] is not None)
    sorted_gcfs = [gid for gid, _ in gcf_counts.most_common()]

    # ── Figure dimensions ───────────────────────────────────────────────────────
    cell_h    = 0.12    # inches per leaf
    tree_w    = 3.0     # tree panel
    strip_w   = 0.18    # width of coupling enzyme strip
    gcf_lbl_w = 0.55    # width of GCF text column
    gap       = 0.05    # gap between panels
    label_w   = 3.2     # genome label area
    left_pad  = 0.2
    right_pad = 1.6     # legend (no GCF color legend needed)

    fig_h = n_leaves * cell_h + 0.5
    fig_w = left_pad + tree_w + gap + strip_w + gap + gcf_lbl_w + gap + label_w + right_pad

    fig = plt.figure(figsize=(fig_w, fig_h))

    f = lambda x: x / fig_w
    heat_bot  = 0.1 / fig_h
    heat_h    = (n_leaves * cell_h) / fig_h
    tree_left = left_pad / fig_w

    ax_tree = fig.add_axes([tree_left,                                        heat_bot, f(tree_w),     heat_h])
    ax_coup = fig.add_axes([tree_left + f(tree_w + gap),                      heat_bot, f(strip_w),    heat_h])
    ax_gcf  = fig.add_axes([tree_left + f(tree_w + gap + strip_w + gap),      heat_bot, f(gcf_lbl_w),  heat_h])
    ax_lbl  = fig.add_axes([tree_left + f(tree_w + gap + strip_w + gap + gcf_lbl_w + gap), heat_bot, f(label_w), heat_h])

    # ── Draw cladogram ──────────────────────────────────────────────────────────
    draw_cladogram(ax_tree, tree.root)
    ax_tree.set_xlim(-0.2, md + 0.5)
    ax_tree.set_ylim(-0.5, n_leaves - 0.5)
    ax_tree.axis('off')

    # ── Draw strips and labels ──────────────────────────────────────────────────
    for plot_idx, rec_id in enumerate(leaf_order):
        meta     = id_to_meta[rec_id]
        gbk_base = meta['gbk_basename']
        gcf_id   = meta['gcf_id']
        genome   = meta['genome'].replace('_', ' ')
        organism = meta['organism']

        # Coupling enzyme class
        cls        = coupling_classes.get(gbk_base, 'Unknown')
        coup_color = COUPLING_COLORS.get(cls, '#aaaaaa')

        y = plot_idx

        # Coupling strip
        ax_coup.barh(y, 1, height=0.85, color=coup_color, left=0)

        # GCF text label
        gcf_text = f'GCF-{gcf_id}' if gcf_id is not None else '—'
        ax_gcf.text(0.5, y, gcf_text, va='center', ha='center',
                    fontsize=5, family='monospace', color='#222222')

        # Genome label: species name + antiSMASH region ID
        label = f'{genome.replace("_", " ")}  {gbk_base}'
        ax_lbl.text(0.02, y, label, va='center', ha='left',
                    fontsize=5.5, family='monospace', color='#222222')

    for ax in (ax_coup, ax_gcf, ax_lbl):
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, n_leaves - 0.5)
        ax.axis('off')

    # ── Strip headers ───────────────────────────────────────────────────────────
    header_y = 1.0 - (0.15 / fig_h)
    for ax, label in [(ax_coup, 'CE'), (ax_gcf, 'GCF')]:
        ax.text(0.5, 1.005, label, transform=ax.transAxes,
                ha='center', va='bottom', fontsize=6, fontweight='bold')

    # ── Legends ────────────────────────────────────────────────────────────────
    coup_patches = [mpatches.Patch(color=COUPLING_COLORS[c], label=c)
                    for c in COUPLING_ORDER if c in coupling_classes.values()]

    legend_x = tree_left + f(tree_w + gap + strip_w + gap + gcf_lbl_w + gap + label_w + 0.15)
    fig.legend(handles=coup_patches, bbox_to_anchor=(legend_x, 0.65),
               loc='upper left', bbox_transform=fig.transFigure,
               fontsize=6.5, title='Coupling enzyme', title_fontsize=7,
               framealpha=0.9)

    fig.suptitle('Phosphonate BGC biosynthetic tree — all BGCs\n'
                 f'NJ tree from BiG-SCAPE pairwise distances  |  '
                 f'{n_leaves} BGCs  |  {len(sorted_gcfs)} GCFs',
                 fontsize=9, y=0.995, va='top')

    # ── Save ───────────────────────────────────────────────────────────────────
    os.makedirs(outdir, exist_ok=True)
    out_png = os.path.join(outdir, 'all_bgcs_biosynthetic_tree.png')
    out_svg = os.path.join(outdir, 'all_bgcs_biosynthetic_tree.svg')
    fig.savefig(out_png, dpi=180, bbox_inches='tight')
    fig.savefig(out_svg,           bbox_inches='tight')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')
    plt.close(fig)


# ─── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='NJ tree of all BGCs from BiG-SCAPE distances',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--db',                  required=True)
    parser.add_argument('--coupling_annotation', required=True)
    parser.add_argument('--outdir',              required=True)
    parser.add_argument('--cutoff', type=float,  default=0.3)
    args = parser.parse_args()

    print('Loading coupling annotations...')
    coupling_classes = load_coupling_classes(args.coupling_annotation)
    print(f'  {len(coupling_classes)} BGC coupling labels loaded')

    print('Loading BGC data from database...')
    conn = sqlite3.connect(args.db)
    bgc_ids, id_to_meta, distances = load_bgc_data(conn, args.cutoff)
    conn.close()
    print(f'  {len(bgc_ids)} BGCs  |  {len(distances)} pairwise distances')

    tree = build_nj_tree(bgc_ids, distances)

    print('Plotting...')
    plot_tree(tree, bgc_ids, id_to_meta, coupling_classes, args.outdir)
    print('Done.')


if __name__ == '__main__':
    main()
