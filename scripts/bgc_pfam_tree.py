#!/usr/bin/env python3
"""
Build a Jaccard-distance NJ phylogenetic tree of BGCs based on Pfam domain content.

Each BGC is represented as a presence/absence vector of Pfam domains. Pairwise
Jaccard distances are computed, then a Neighbor-Joining tree is constructed.

Usage:
    python scripts/bgc_pfam_tree.py \
        --db results/bigscape_results/Pantoea/Pantoea.db \
        --bgc_type phosphonate \
        --outdir results/bgc_trees/Pantoea

Output files:
    {bgc_type}_pfam_tree.nwk       - Newick tree (load in iTOL, FigTree, etc.)
    {bgc_type}_jaccard_distances.tsv - Full pairwise distance matrix
    {bgc_type}_domain_matrix.tsv   - BGC x Pfam presence/absence matrix
    {bgc_type}_metadata.json       - BGC labels, products, GCF family assignments
"""

import argparse
import sqlite3
import os
import sys
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent))
from utils.bgc_labels import label_from_path, make_labels_unique
from utils.tree_building import build_nj_tree

# NJ trees for large datasets can be deeply nested; raise the limit
sys.setrecursionlimit(10000)


# ─── Database queries ────────────────────────────────────────────────────────

BGC_QUERY = """
    SELECT
        br.id       AS bgc_id,
        br.product,
        g.id        AS gbk_id,
        g.path      AS gbk_path,
        g.organism
    FROM bgc_record br
    JOIN gbk g ON br.gbk_id = g.id
    WHERE LOWER(br.product) LIKE ?
    ORDER BY g.path
"""

BGC_QUERY_FAMILY = """
    SELECT
        br.id       AS bgc_id,
        br.product,
        g.id        AS gbk_id,
        g.path      AS gbk_path,
        g.organism
    FROM bgc_record br
    JOIN gbk g ON br.gbk_id = g.id
    JOIN bgc_record_family brf ON br.id = brf.record_id
    WHERE LOWER(br.product) LIKE ?
      AND brf.family_id = ?
    ORDER BY g.path
"""

DOMAIN_QUERY = """
    SELECT DISTINCT h.accession
    FROM cds c
    JOIN scanned_cds sc ON sc.cds_id = c.id
    JOIN hsp h ON h.cds_id = c.id
    WHERE c.gbk_id = ?
      AND h.accession != ''
      AND h.bit_score >= 20
"""

FAMILY_QUERY = """
    SELECT brf.family_id, f.cutoff
    FROM bgc_record_family brf
    JOIN family f ON brf.family_id = f.id
    WHERE brf.record_id = ?
    ORDER BY f.cutoff
"""


# ─── Core logic ──────────────────────────────────────────────────────────────

def load_bgc_domains(db_path, bgc_type_filter, family_id=None):
    """
    Query BiG-SCAPE DB for all BGCs matching bgc_type_filter and their
    Pfam domain sets. Optionally restrict to a single GCF family.

    Returns:
        bgc_labels   : list of str, one per BGC
        bgc_domains  : list of set, Pfam accessions per BGC
        bgc_metadata : list of dict with product, organism, gcf info
    """
    bgc_labels   = []
    bgc_domains  = []
    bgc_metadata = []
    skipped = 0

    with sqlite3.connect(db_path) as conn:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()

        if family_id is not None:
            cur.execute(BGC_QUERY_FAMILY, (f'%{bgc_type_filter.lower()}%', family_id))
            print(f"  Filtering to GCF family {family_id}")
        else:
            cur.execute(BGC_QUERY, (f'%{bgc_type_filter.lower()}%',))
        bgc_rows = cur.fetchall()

        if not bgc_rows:
            raise ValueError(f"No BGCs found matching '{bgc_type_filter}' in {db_path}")

        print(f"  Found {len(bgc_rows)} BGCs matching '{bgc_type_filter}'")

        for row in bgc_rows:
            # Get Pfam domains for this BGC's GBK
            cur.execute(DOMAIN_QUERY, (row['gbk_id'],))
            domains = {r[0] for r in cur.fetchall()}

            if not domains:
                skipped += 1
                continue

            # Get GCF family assignment(s) for metadata
            cur.execute(FAMILY_QUERY, (row['bgc_id'],))
            families = [{'family_id': r['family_id'], 'cutoff': r['cutoff']}
                        for r in cur.fetchall()]

            label = label_from_path(row['gbk_path'])
            bgc_labels.append(label)
            bgc_domains.append(domains)
            bgc_metadata.append({
                'label':     label,
                'product':   row['product'],
                'organism':  row['organism'] or os.path.basename(os.path.dirname(row['gbk_path'])),
                'gbk_path':  row['gbk_path'],
                'n_domains': len(domains),
                'families':  families,
            })

    if skipped:
        print(f"  Skipped {skipped} BGCs with no Pfam domain hits")

    bgc_labels = make_labels_unique(bgc_labels)
    for i, meta in enumerate(bgc_metadata):
        meta['label'] = bgc_labels[i]

    return bgc_labels, bgc_domains, bgc_metadata


def build_domain_matrix(bgc_labels, bgc_domains):
    """
    Build a binary BGC × Pfam presence/absence matrix.
    Returns: all_domains (sorted list), numpy array (n_bgc x n_domains)
    """
    all_domains = sorted(set().union(*bgc_domains))
    domain_index = {d: i for i, d in enumerate(all_domains)}
    n = len(bgc_labels)
    m = len(all_domains)

    matrix = np.zeros((n, m), dtype=np.uint8)
    for i, domains in enumerate(bgc_domains):
        for d in domains:
            matrix[i, domain_index[d]] = 1

    return all_domains, matrix


def build_jaccard_matrix(matrix):
    """
    Compute pairwise Jaccard distances between BGC vectors.
    Jaccard distance = 1 - (|A ∩ B| / |A ∪ B|)
    Uses scipy cdist for vectorised computation.
    Returns: numpy (n x n) distance matrix
    """
    dist = cdist(matrix.astype(float), matrix.astype(float), metric='jaccard')
    np.fill_diagonal(dist, 0.0)
    return dist


# ─── Output ───────────────────────────────────────────────────────────────────

def save_outputs(outdir, bgc_type, bgc_labels, all_domains, domain_matrix,
                 dist_matrix, tree, bgc_metadata):
    os.makedirs(outdir, exist_ok=True)

    # Newick tree
    newick_path = os.path.join(outdir, f"{bgc_type}_pfam_tree.nwk")
    Phylo.write(tree, newick_path, 'newick')
    print(f"  Tree:            {newick_path}")

    # Distance matrix TSV
    dist_path = os.path.join(outdir, f"{bgc_type}_jaccard_distances.tsv")
    with open(dist_path, 'w') as f:
        f.write('\t' + '\t'.join(bgc_labels) + '\n')
        for i, label in enumerate(bgc_labels):
            row = '\t'.join(f"{dist_matrix[i][j]:.4f}" for j in range(len(bgc_labels)))
            f.write(f"{label}\t{row}\n")
    print(f"  Distance matrix: {dist_path}")

    # Domain presence/absence matrix TSV
    domain_path = os.path.join(outdir, f"{bgc_type}_domain_matrix.tsv")
    with open(domain_path, 'w') as f:
        f.write('\t' + '\t'.join(all_domains) + '\n')
        for i, label in enumerate(bgc_labels):
            row = '\t'.join(str(domain_matrix[i][j]) for j in range(len(all_domains)))
            f.write(f"{label}\t{row}\n")
    print(f"  Domain matrix:   {domain_path}")

    # Metadata JSON
    meta_path = os.path.join(outdir, f"{bgc_type}_metadata.json")
    with open(meta_path, 'w') as f:
        json.dump(bgc_metadata, f, indent=2)
    print(f"  Metadata:        {meta_path}")

    # Quick ASCII tree preview in terminal
    print("\n--- Tree preview (ASCII) ---")
    buf = StringIO()
    Phylo.draw_ascii(tree, file=buf)
    lines = buf.getvalue().splitlines()
    # Print first 40 lines to avoid flooding the terminal
    for line in lines[:40]:
        print(line)
    if len(lines) > 40:
        print(f"  ... ({len(lines) - 40} more lines — open the .nwk file for the full tree)")


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build a Jaccard-distance NJ tree of BGCs from Pfam domain content',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--db',       required=True,
                        help='Path to BiG-SCAPE SQLite database (.db)')
    parser.add_argument('--bgc_type', default='phosphonate',
                        help='BGC product type filter, case-insensitive substring match '
                             '(default: phosphonate)')
    parser.add_argument('--outdir',   required=True,
                        help='Output directory')
    parser.add_argument('--family_id', type=int, default=None,
                        help='Restrict tree to a single GCF family ID (e.g. --family_id 2)')
    parser.add_argument('--metadata_only', action='store_true',
                        help='Only write metadata JSON; skip tree/matrix/distance outputs '
                             '(much faster for large datasets)')
    args = parser.parse_args()

    print(f"Loading BGC Pfam domains from: {args.db}")
    bgc_labels, bgc_domains, bgc_metadata = load_bgc_domains(
        args.db, args.bgc_type, family_id=args.family_id)

    if args.metadata_only:
        os.makedirs(args.outdir, exist_ok=True)
        meta_path = os.path.join(args.outdir, f"{args.bgc_type}_metadata.json")
        with open(meta_path, 'w') as f:
            json.dump(bgc_metadata, f, indent=2)
        print(f"  Metadata: {meta_path}")
        print("Done (metadata only).")
        return

    if len(bgc_labels) < 3:
        raise ValueError(f"Need at least 3 BGCs to build a tree, found {len(bgc_labels)}")

    print(f"Building domain matrix ({len(bgc_labels)} BGCs)...")
    all_domains, domain_matrix = build_domain_matrix(bgc_labels, bgc_domains)
    print(f"  Unique Pfam domains: {len(all_domains)}")

    print("Computing Jaccard distances...")
    dist_matrix = build_jaccard_matrix(domain_matrix)

    print("Constructing NJ tree...")
    tree = build_nj_tree(bgc_labels, dist_matrix)

    print("Saving outputs...")
    save_outputs(args.outdir, args.bgc_type, bgc_labels, all_domains,
                 domain_matrix, dist_matrix, tree, bgc_metadata)

    print("\nDone. Load the .nwk file in iTOL (https://itol.embl.de) or FigTree for visualization.")


if __name__ == '__main__':
    main()
