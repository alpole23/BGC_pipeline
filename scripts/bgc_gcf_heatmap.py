#!/usr/bin/env python3
"""
GCF × species presence/absence heatmap for phosphonate BGCs.

Rows (clustered by profile similarity):
  - One per GCF, labelled "GCF-X (n=total)"
  - One "Singletons" row per coupling enzyme class that has singletons
  - Ordered by hierarchical clustering of species co-occurrence profiles
  - Row dendrogram shown on the left

Columns:
  - One per species, ordered by GTDB-Tk phylogenetic tree (if provided),
    otherwise sorted alphabetically by genus then epithet
  - Column phylogenetic tree shown above the heatmap

Cell color:
  - Coupling enzyme class color (lightened) if present, white if absent
  - Cell text shows count

Left strip:
  - Coupling enzyme class color per row

Usage (with phylogenetic tree):
    python scripts/bgc_gcf_heatmap.py \\
        --metadata            results/bgc_trees/Pantoea/phosphonate_metadata.json \\
        --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \\
        --gtdbtk_tree         results/gtdbtk_results/Pantoea/gtdbtk_output/classify/gtdbtk.bac120.classify.tree.1.tree \\
        --gtdbtk_summary      results/gtdbtk_results/Pantoea/gtdbtk_output/gtdbtk.bac120.summary.tsv \\
        --outdir              results/bgc_trees/Pantoea
"""

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import (COUPLING_COLORS, COUPLING_ORDER as CLASS_ORDER,
                              LEGACY_CLASS_NAMES as LEGACY)


# ─── GCF biosynthetic tree helpers ────────────────────────────────────────────

def load_gcf_tree(path):
    """Read the NJ GCF biosynthetic tree (Newick) from bgc_gcf_tree.py."""
    return Phylo.read(path, 'newick')


def order_rows_by_tree(rows, tree):
    """
    Reorder GCF rows to match the NJ tree leaf traversal order.
    Singleton rows are inserted after the last GCF of their coupling class.
    GCFs absent from the tree are appended at the end.
    """
    leaf_ids = []
    for clade in tree.get_terminals():
        name = clade.name.strip("'") if clade.name else ''
        try:
            leaf_ids.append(int(name.replace('GCF-', '')))
        except ValueError:
            pass

    gcf_rows      = {row['gcf_id']: row for row in rows if row['type'] == 'gcf'}
    singleton_rows = {row['cls']: row   for row in rows if row['type'] == 'singleton'}

    ordered   = []
    last_cls  = None
    used_sing = set()

    for gcf_id in leaf_ids:
        if gcf_id not in gcf_rows:
            continue
        row = gcf_rows[gcf_id]
        # Insert singleton for previous class when class boundary is crossed
        if last_cls is not None and row['cls'] != last_cls:
            if last_cls in singleton_rows and last_cls not in used_sing:
                ordered.append(singleton_rows[last_cls])
                used_sing.add(last_cls)
        ordered.append(row)
        last_cls = row['cls']

    # Singleton for the last class in the tree
    if last_cls and last_cls in singleton_rows and last_cls not in used_sing:
        ordered.append(singleton_rows[last_cls])
        used_sing.add(last_cls)

    # Append any GCFs not in the tree, then remaining singletons
    tree_gcf_ids = set(leaf_ids)
    for row in rows:
        if row['type'] == 'gcf' and row['gcf_id'] not in tree_gcf_ids:
            ordered.append(row)
    for cls, row in singleton_rows.items():
        if cls not in used_sing:
            ordered.append(row)

    return ordered


def assign_row_layout(clade, row_positions, depth=0):
    """
    Assign ._depth and ._x (y-position in plot) to every node of the GCF tree.
    row_positions: {gcf_id: row_index} for terminal nodes.
    Internal nodes get ._x as the mean of their children's ._x values.
    """
    clade._depth = depth
    if clade.is_terminal():
        name   = clade.name.strip("'") if clade.name else ''
        try:
            gcf_id = int(name.replace('GCF-', ''))
        except ValueError:
            gcf_id = -1
        clade._x = row_positions.get(gcf_id, 0)
        return
    for child in clade.clades:
        assign_row_layout(child, row_positions, depth + 1)
    clade._x = sum(c._x for c in clade.clades) / len(clade.clades)


def _draw_row_cladogram(ax, clade, color='#333333', lw=1.2):
    """
    Draw the GCF cladogram on a vertical axis.
    Root on the LEFT (low x = low depth), leaves on the RIGHT (high depth),
    y-axis matches the heatmap row positions.
    """
    if clade.is_terminal():
        return
    x_node    = clade._depth
    child_ys  = [c._x for c in clade.clades]
    # Vertical bar connecting children at this node's depth
    ax.plot([x_node, x_node], [min(child_ys), max(child_ys)],
            color=color, lw=lw, solid_capstyle='round')
    for child in clade.clades:
        # Horizontal line from node to child
        ax.plot([x_node, child._depth], [child._x, child._x],
                color=color, lw=lw, solid_capstyle='round')
        _draw_row_cladogram(ax, child, color, lw)


# ─── Data loading ─────────────────────────────────────────────────────────────

def load_coupling_classes(path):
    classes = {}
    in_data = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == 'DATA':
                in_data = True
                continue
            if in_data and line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    cls = LEGACY.get(parts[2], parts[2])
                    classes[parts[0]] = cls
    return classes


def species_from_organism(organism):
    """Extract 'Genus species' from organism string (first two words)."""
    parts = organism.strip().split()
    if len(parts) >= 2:
        return f'{parts[0]} {parts[1]}'
    return organism.strip() or 'Unknown'


# ─── Build data matrix ────────────────────────────────────────────────────────

def build_matrix(metadata_path, coupling_classes):
    with open(metadata_path) as f:
        meta = json.load(f)

    # Per-BGC records — only region-level records (those with a GCF assignment at 0.3)
    # Sub-records (cand_cluster, protocluster) have no GCF and are excluded.
    records = []
    for r in meta:
        gcf = next(
            (f['family_id'] for f in r.get('families', []) if f['cutoff'] == 0.3),
            None,
        )
        if gcf is None:
            continue  # skip cand_cluster / protocluster sub-records
        lbl     = r['label']
        species = species_from_organism(r.get('organism', ''))
        cls     = coupling_classes.get(lbl, 'Unknown')
        records.append({'label': lbl, 'species': species, 'gcf': gcf, 'cls': cls,
                        'gbk_path': r.get('gbk_path', '')})

    # GCF sizes (total BGCs regardless of species)
    gcf_total = Counter(rec['gcf'] for rec in records)

    # Assign coupling class to each GCF by majority vote
    gcf_cls_votes = defaultdict(Counter)
    for rec in records:
        gcf_cls_votes[rec['gcf']][rec['cls']] += 1
    gcf_class = {gcf: votes.most_common(1)[0][0] for gcf, votes in gcf_cls_votes.items()}

    # Single-member GCFs are treated as "singletons" (isolated BGCs)
    singleton_gcf_ids = {gid for gid, cnt in gcf_total.items() if cnt == 1}

    # Species list (sorted: genus then epithet)
    all_species = sorted(set(rec['species'] for rec in records),
                         key=lambda s: (s.split()[0], s))

    # GCF × species count matrix (multi-member GCFs only)
    gcf_species = defaultdict(Counter)
    for rec in records:
        if rec['gcf'] not in singleton_gcf_ids:
            gcf_species[rec['gcf']][rec['species']] += 1

    # Singleton × class × species count matrix (single-member GCFs)
    singleton_species = defaultdict(Counter)
    for rec in records:
        if rec['gcf'] in singleton_gcf_ids:
            singleton_species[rec['cls']][rec['species']] += 1

    n_singletons_total = sum(sum(v.values()) for v in singleton_species.values())

    # Remove singleton GCFs from gcf_total and gcf_class
    for gid in singleton_gcf_ids:
        del gcf_total[gid]
        del gcf_class[gid]

    return (records, gcf_class, gcf_total, gcf_species,
            singleton_species, all_species, n_singletons_total)


# ─── Build ordered row list ───────────────────────────────────────────────────

def build_rows(gcf_class, gcf_total, gcf_species, singleton_species):
    """
    Returns list of row dicts, ordered by coupling class then GCF size.
    Each dict: {label, cls, type ('gcf'|'singleton'), gcf_id (or None), total}
    """
    rows = []
    for cls in CLASS_ORDER:
        # GCFs in this class, sorted by total size descending
        class_gcfs = sorted(
            [(gid, gcf_total[gid]) for gid, c in gcf_class.items() if c == cls],
            key=lambda x: -x[1],
        )
        for gid, total in class_gcfs:
            rows.append({
                'label':  f'GCF-{gid}  (n={total})',
                'cls':    cls,
                'type':   'gcf',
                'gcf_id': gid,
                'total':  total,
            })
        # Singleton row for this class (if any exist)
        if singleton_species.get(cls):
            n = sum(singleton_species[cls].values())
            rows.append({
                'label':  f'Singletons  (n={n})',
                'cls':    cls,
                'type':   'singleton',
                'gcf_id': None,
                'total':  n,
            })
    return rows


# ─── Row clustering ───────────────────────────────────────────────────────────

def compute_row_linkage(rows, species_list, gcf_species, singleton_species):
    """
    Hierarchical clustering of rows by species presence/absence (Jaccard).
    Returns (linkage_matrix, reordered_row_indices).
    """
    n_rows = len(rows)
    n_cols = len(species_list)
    sp_idx = {sp: i for i, sp in enumerate(species_list)}

    matrix = np.zeros((n_rows, n_cols), dtype=float)
    for ri, row in enumerate(rows):
        if row['type'] == 'gcf':
            cell_data = gcf_species[row['gcf_id']]
        else:
            cell_data = singleton_species[row['cls']]
        for sp, cnt in cell_data.items():
            ci = sp_idx.get(sp)
            if ci is not None:
                matrix[ri, ci] = 1.0

    if n_rows < 2:
        return None, list(range(n_rows))

    dist = pdist(matrix, metric='jaccard')
    dist = np.nan_to_num(dist, nan=1.0)
    Z = linkage(dist, method='complete')
    order = list(leaves_list(Z))
    return Z, order


def draw_row_dendrogram(ax, Z, n_rows):
    """
    Draw row dendrogram in ax. Leaves on right (adjacent to strip), root on left.
    Y-axis matches imshow convention: row 0 at top.
    """
    dend = dendrogram(Z, no_plot=True, orientation='left')
    # scipy leaf positions: 5, 15, 25, ... (5 + 10*i)
    # Map to row positions 0..n-1: row_pos = (y - 5) / 10
    max_dist = max(d for dc in dend['dcoord'] for d in dc)
    for ic, dc in zip(dend['icoord'], dend['dcoord']):
        y_mapped = [(y - 5) / 10 for y in ic]
        ax.plot(dc, y_mapped, color='#333333', lw=0.9, solid_capstyle='butt')
    # Match imshow y-axis: row 0 near top, row n-1 near bottom
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.invert_xaxis()   # leaves on right, root on left
    ax.set_xlim(max_dist * 1.05, 0)
    ax.axis('off')


# ─── Phylogenetic tree ────────────────────────────────────────────────────────

def _keep_clade(clade, targets):
    """Prune in-place. Returns True if clade has ≥1 target descendant."""
    if clade.is_terminal():
        return clade.name in targets
    clade.clades = [c for c in clade.clades if _keep_clade(c, targets)]
    return len(clade.clades) > 0


def _assign_layout(clade, counter, depth=0):
    """Assign _x (leaf traversal position) and _depth to every node."""
    clade._depth = depth
    if clade.is_terminal():
        clade._x = counter[0]
        counter[0] += 1
        return [clade]
    leaves = []
    for child in clade.clades:
        leaves.extend(_assign_layout(child, counter, depth + 1))
    clade._x = (clade.clades[0]._x + clade.clades[-1]._x) / 2.0
    return leaves


def _max_depth(clade):
    if clade.is_terminal():
        return clade._depth
    return max(_max_depth(c) for c in clade.clades)


def _draw_cladogram(ax, clade, max_depth, color='#333333', lw=1.2):
    """Draw cladogram: root at top (y=max_depth), leaves at bottom (y=0)."""
    if clade.is_terminal():
        return
    y_node = max_depth - clade._depth
    child_xs = [c._x for c in clade.clades]
    ax.plot([min(child_xs), max(child_xs)], [y_node, y_node],
            color=color, lw=lw, solid_capstyle='round')
    for child in clade.clades:
        y_child = max_depth - child._depth
        ax.plot([child._x, child._x], [y_node, y_child],
                color=color, lw=lw, solid_capstyle='round')
        _draw_cladogram(ax, child, max_depth, color, lw)


def load_and_prune_phylo_tree(tree_path, summary_path, records):
    """
    Prune GTDB-Tk tree to one representative genome per heatmap species.
    Returns (pruned_tree, ordered_species_list) or (None, None) on failure.
    """
    # genome → heatmap organism name (from metadata)
    org_genomes = defaultdict(list)
    for rec in records:
        gbk = rec.get('gbk_path', '')
        if not gbk:
            continue
        genome = os.path.basename(os.path.dirname(gbk))
        org_genomes[rec['species']].append(genome)

    # Which genomes are in GTDB-Tk?
    gtdbtk_set = set()
    with open(summary_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            gtdbtk_set.add(row['user_genome'].replace('usr_', '', 1))

    # Pick first GTDB-Tk genome per organism group as tree representative
    # target_leaves: 'usr_genome' → organism_name
    target_leaves = {}
    for org, genomes in org_genomes.items():
        for g in sorted(set(genomes)):
            if g in gtdbtk_set:
                target_leaves[f'usr_{g}'] = org
                break

    if not target_leaves:
        print('  Warning: no GTDB-Tk representatives found — skipping phylo tree')
        return None, None

    print(f'  Pruning GTDB-Tk tree to {len(target_leaves)} representative genomes...')
    tree = Phylo.read(tree_path, 'newick')
    _keep_clade(tree.root, set(target_leaves.keys()))

    # Relabel leaves with organism names
    for clade in tree.get_terminals():
        if clade.name in target_leaves:
            clade.name = target_leaves[clade.name]

    # Compute layout and get leaf order
    _assign_layout(tree.root, [0])
    ordered = [c.name for c in tree.get_terminals()]

    placed   = set(ordered)
    missing  = [sp for sp in org_genomes if sp not in placed]
    if missing:
        print(f'  Warning: {len(missing)} species not placed in tree: {missing}')

    return tree, ordered


# ─── Color helpers ────────────────────────────────────────────────────────────

def hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))


def lighten(hex_color, amount=0.25):
    r, g, b = hex_to_rgb(hex_color)
    return (r + (1-r)*amount, g + (1-g)*amount, b + (1-b)*amount)


# ─── Plot ─────────────────────────────────────────────────────────────────────

def plot_heatmap(rows, species_list, gcf_species, singleton_species, outdir,
                 phylo_tree=None, phylo_species_order=None, gcf_tree=None):

    # ── Row order ───────────────────────────────────────────────────────────────
    if gcf_tree is not None:
        rows = order_rows_by_tree(rows, gcf_tree)
    # else: rows keep the coupling-class order from build_rows()

    # ── Column order ────────────────────────────────────────────────────────────
    if phylo_tree is not None and phylo_species_order:
        # Species order from phylogenetic tree, unplaced species appended at end
        placed = set(phylo_species_order)
        extra  = [sp for sp in species_list if sp not in placed]
        species_list = phylo_species_order + extra
    # else: keep the alphabetical species_list from build_matrix

    n_rows = len(rows)
    n_cols = len(species_list)

    # ── Build image and count arrays ────────────────────────────────────────────
    img    = np.ones((n_rows, n_cols, 4))   # RGBA, white by default
    counts = np.zeros((n_rows, n_cols), dtype=int)
    sp_idx = {sp: i for i, sp in enumerate(species_list)}

    for ri, row in enumerate(rows):
        cls   = row['cls']
        light = lighten(COUPLING_COLORS.get(cls, '#aaaaaa'), amount=0.35)
        cell_data = (gcf_species[row['gcf_id']] if row['type'] == 'gcf'
                     else singleton_species[cls])
        for sp, cnt in cell_data.items():
            ci = sp_idx.get(sp)
            if ci is not None:
                counts[ri, ci] = cnt
                img[ri, ci, :3] = light
                img[ri, ci,  3] = 1.0

    # ── Figure dimensions ───────────────────────────────────────────────────────
    dendro_w  = 1.5 if gcf_tree is not None else 0.0  # GCF biosynthetic tree
    strip_w   = 0.25    # coupling class color strip
    cell_w    = 1.0     # per species column
    cell_h    = 0.45    # per row
    label_w   = 2.8     # row label area (right of strip)
    bot_pad   = 0.4
    right_pad = 2.5     # space for legend

    has_tree  = phylo_tree is not None
    tree_h    = 2.5 if has_tree else 0.0
    top_pad   = 2.5 if has_tree else 2.5   # space for rotated col labels (+ tree if present)

    fig_w = dendro_w + strip_w + label_w + n_cols * cell_w + right_pad
    fig_h = top_pad + tree_h + n_rows * cell_h + bot_pad

    fig = plt.figure(figsize=(fig_w, fig_h))

    # ── Compute figure-fraction positions ───────────────────────────────────────
    f_dendro_w = dendro_w  / fig_w
    f_strip_w  = strip_w   / fig_w
    f_label_w  = label_w   / fig_w
    f_heat_w   = (n_cols * cell_w) / fig_w
    f_heat_h   = (n_rows * cell_h) / fig_h
    f_tree_h   = tree_h / fig_h
    f_bot      = bot_pad / fig_h
    f_heat_bot = f_bot
    f_tree_bot = f_bot + f_heat_h

    # Strip is placed immediately left of heatmap; labels flow left into the gap
    f_strip_left = f_dendro_w + f_label_w
    f_heat_left  = f_dendro_w + f_label_w + f_strip_w

    # ── Create axes ─────────────────────────────────────────────────────────────
    ax_row_dendro = fig.add_axes([0,            f_heat_bot, f_dendro_w, f_heat_h])
    ax_strip      = fig.add_axes([f_strip_left, f_heat_bot, f_strip_w,  f_heat_h])
    ax_heat       = fig.add_axes([f_heat_left,  f_heat_bot, f_heat_w,   f_heat_h])

    if has_tree:
        ax_col_tree = fig.add_axes([f_heat_left, f_tree_bot, f_heat_w, f_tree_h])

    # ── Row GCF biosynthetic tree ────────────────────────────────────────────────
    if gcf_tree is not None and dendro_w > 0:
        gcf_row_positions = {row['gcf_id']: ri
                             for ri, row in enumerate(rows) if row['type'] == 'gcf'}
        assign_row_layout(gcf_tree.root, gcf_row_positions)
        md = _max_depth(gcf_tree.root)
        _draw_row_cladogram(ax_row_dendro, gcf_tree.root)
        ax_row_dendro.set_xlim(-0.1, md + 0.3)
        ax_row_dendro.set_ylim(n_rows - 0.5, -0.5)
    ax_row_dendro.axis('off')

    # ── Class color strip + row labels ──────────────────────────────────────────
    strip_img = np.ones((n_rows, 1, 4))
    for ri, row in enumerate(rows):
        c = hex_to_rgb(COUPLING_COLORS.get(row['cls'], '#aaaaaa'))
        strip_img[ri, 0, :3] = c
        strip_img[ri, 0,  3] = 1.0

    ax_strip.imshow(strip_img, aspect='auto', interpolation='nearest')
    ax_strip.set_xticks([])
    ax_strip.set_yticks(range(n_rows))
    ax_strip.set_yticklabels([row['label'] for row in rows],
                              fontsize=8.5, va='center')
    ax_strip.tick_params(axis='y', length=0, pad=6)

    # Separator lines between coupling classes (on both axes)
    prev_cls = None
    for ri, row in enumerate(rows):
        if prev_cls is not None and row['cls'] != prev_cls:
            for ax in (ax_strip, ax_heat):
                ax.axhline(ri - 0.5, color='white', lw=2, zorder=3)
        prev_cls = row['cls']

    # ── Main heatmap ─────────────────────────────────────────────────────────────
    ax_heat.imshow(img, aspect='auto', interpolation='nearest')

    for x in range(n_cols + 1):
        ax_heat.axvline(x - 0.5, color='#dddddd', lw=0.5, zorder=2)
    for y in range(n_rows + 1):
        ax_heat.axhline(y - 0.5, color='#dddddd', lw=0.5, zorder=2)

    for ri in range(n_rows):
        for ci in range(n_cols):
            cnt = counts[ri, ci]
            if cnt > 0:
                ax_heat.text(ci, ri, str(cnt),
                             ha='center', va='center', fontsize=7,
                             color='#222222', fontweight='bold')

    ax_heat.set_yticks([])

    if has_tree:
        # Column labels go on the tree axes (top)
        ax_heat.set_xticks([])
    else:
        ax_heat.set_xticks(range(n_cols))
        ax_heat.set_xticklabels(species_list, rotation=45, ha='right',
                                 fontsize=9, style='italic')
        ax_heat.xaxis.set_ticks_position('top')
        ax_heat.xaxis.set_label_position('top')
        ax_heat.tick_params(axis='x', length=0, pad=3)

    # ── Column phylogenetic tree ─────────────────────────────────────────────────
    if has_tree:
        md = _max_depth(phylo_tree.root)
        _draw_cladogram(ax_col_tree, phylo_tree.root, md)

        ax_col_tree.set_xlim(-0.5, n_cols - 0.5)
        ax_col_tree.set_ylim(-0.15, md + 0.3)
        ax_col_tree.set_xticks(range(n_cols))
        ax_col_tree.set_xticklabels(species_list, rotation=45, ha='left',
                                     va='bottom', fontsize=9, style='italic')
        ax_col_tree.xaxis.set_ticks_position('top')
        ax_col_tree.xaxis.set_label_position('top')
        ax_col_tree.tick_params(axis='x', length=0, pad=3)
        ax_col_tree.set_yticks([])
        for spine in ax_col_tree.spines.values():
            spine.set_visible(False)

    # ── Legend ──────────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=COUPLING_COLORS[cls], label=cls)
        for cls in CLASS_ORDER
        if cls in {row['cls'] for row in rows}
    ]
    # Place legend in the right margin (explicit figure-fraction coordinates)
    legend_x = f_heat_left + f_heat_w + 0.15 / fig_w   # just right of heatmap
    legend_y = 0.5
    fig.legend(handles=legend_patches,
               bbox_to_anchor=(legend_x, legend_y),
               loc='center left',
               bbox_transform=fig.transFigure,
               fontsize=8, title='Coupling enzyme class',
               title_fontsize=8.5, framealpha=0.9)

    # ── Title ───────────────────────────────────────────────────────────────────
    n_gcfs      = sum(1 for r in rows if r['type'] == 'gcf')
    n_sing_rows = sum(1 for r in rows if r['type'] == 'singleton')
    n_sing_bgcs = sum(r['total'] for r in rows if r['type'] == 'singleton')
    fig.suptitle(
        f'Phosphonate BGC diversity — GCF × Species\n'
        f'{n_gcfs} GCFs  |  {n_sing_bgcs} singletons ({n_sing_rows} classes)  |  {n_cols} species',
        fontsize=11, y=0.99, va='top',
    )

    # ── Save ────────────────────────────────────────────────────────────────────
    os.makedirs(outdir, exist_ok=True)
    out_png = os.path.join(outdir, 'gcf_species_heatmap.png')
    out_svg = os.path.join(outdir, 'gcf_species_heatmap.svg')
    fig.savefig(out_png, dpi=180, bbox_inches='tight')
    fig.savefig(out_svg,           bbox_inches='tight')
    print(f'Saved: {out_png}')
    print(f'Saved: {out_svg}')
    plt.close(fig)


# ─── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='GCF × species heatmap for phosphonate BGCs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--metadata',            required=True)
    parser.add_argument('--coupling_annotation', required=True)
    parser.add_argument('--outdir',              required=True)
    parser.add_argument('--gtdbtk_tree',   default=None,
                        help='GTDB-Tk classify tree (Newick) for column ordering')
    parser.add_argument('--gtdbtk_summary', default=None,
                        help='GTDB-Tk summary TSV for genome→species mapping')
    parser.add_argument('--gcf_tree', default=None,
                        help='GCF biosynthetic NJ tree (Newick) from bgc_gcf_tree.py '
                             'for row ordering')
    args = parser.parse_args()

    print('Loading data...')
    coupling_classes = load_coupling_classes(args.coupling_annotation)

    (records, gcf_class, gcf_total, gcf_species,
     singleton_species, all_species, n_singletons) = build_matrix(
        args.metadata, coupling_classes)

    print(f'  {len(all_species)} species  |  {len(gcf_total)} GCFs  |  {n_singletons} singletons')
    print('  Coupling class → GCF mapping:')
    for gid in sorted(gcf_class):
        print(f'    GCF-{gid}: {gcf_class[gid]} (n={gcf_total[gid]})')

    rows = build_rows(gcf_class, gcf_total, gcf_species, singleton_species)
    print(f'  {len(rows)} rows ({sum(1 for r in rows if r["type"]=="gcf")} GCFs + '
          f'{sum(1 for r in rows if r["type"]=="singleton")} singleton rows)')

    # Phylogenetic tree for column ordering
    phylo_tree = phylo_order = None
    if args.gtdbtk_tree and args.gtdbtk_summary:
        print('Loading phylogenetic tree...')
        phylo_tree, phylo_order = load_and_prune_phylo_tree(
            args.gtdbtk_tree, args.gtdbtk_summary, records)

    # GCF biosynthetic tree for row ordering
    gcf_tree = None
    if args.gcf_tree:
        print('Loading GCF biosynthetic tree...')
        gcf_tree = load_gcf_tree(args.gcf_tree)

    print('Plotting...')
    plot_heatmap(rows, all_species, gcf_species, singleton_species, args.outdir,
                 phylo_tree=phylo_tree, phylo_species_order=phylo_order,
                 gcf_tree=gcf_tree)
    print('Done.')


if __name__ == '__main__':
    main()
