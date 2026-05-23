#!/usr/bin/env python3
"""
Build a gene-order (synteny) NJ tree of BGCs.

Each BGC is represented as an ORDERED sequence of Pfam domains — one per
annotated gene, sorted by genomic position. Pairwise distances are computed
using normalized LCS (Longest Common Subsequence), which captures shared
pathway architecture while tolerating gene insertions, deletions, and
duplications.

    distance(A, B) = 1 - LCS(A, B) / max(|A|, |B|)

Example:
    BGCa: pepM - monooxygenase - dehydratase - ATPgrasp - reductase
    BGCb: pepM - monooxygenase - dehydratase - reductase
    LCS  = [pepM, monooxygenase, dehydratase, reductase]  →  4 shared
    dist = 1 - 4/5 = 0.20  (close)

    BGCa_GCFa: pepM - monooxygenase - dehydratase - ATPgrasp - reductase
    BGCa_GCFb: pepM - decarboxylase - aminotransferase - ATPgrasp
    LCS  = [pepM, ATPgrasp]  →  2 shared
    dist = 1 - 2/5 = 0.60  (distant — diverge after pepM)

Usage:
    # Single GCF
    python scripts/bgc_synteny_tree.py \\
        --db results/bigscape_results/Pantoea/Pantoea.db \\
        --bgc_type phosphonate \\
        --family_id 2 \\
        --outdir results/bgc_trees/Pantoea/GCF2_synteny

    # All phosphonate BGCs (cross-GCF comparison)
    python scripts/bgc_synteny_tree.py \\
        --db results/bigscape_results/Pantoea/Pantoea.db \\
        --bgc_type phosphonate \\
        --outdir results/bgc_trees/Pantoea/all_synteny

Outputs:
    {bgc_type}_synteny_tree.nwk        Newick tree
    {bgc_type}_domain_sequences.tsv    Ordered domain sequence per BGC
    {bgc_type}_lcs_distances.tsv       Pairwise distance matrix
    {bgc_type}_metadata.json           BGC metadata
"""

import argparse
import sqlite3
import os
import sys
import json
from pathlib import Path

import numpy as np
from Bio import Phylo
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import DOMAIN_NAMES
from utils.bgc_labels import label_from_path, make_labels_unique
from utils.tree_building import build_nj_tree

sys.setrecursionlimit(10000)


def domain_name(accession):
    """Return a readable name for a Pfam accession, stripping version suffix."""
    base = accession.split('.')[0]
    return DOMAIN_NAMES.get(base, base)


# ─── Database queries ─────────────────────────────────────────────────────────

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

# Returns all domain hits per CDS, sorted by genomic position then bit score.
# We group in Python to take the best hit per CDS.
DOMAIN_SEQ_QUERY = """
    SELECT
        c.id        AS cds_id,
        c.nt_start,
        c.strand,
        h.accession,
        h.bit_score
    FROM cds c
    JOIN scanned_cds sc ON sc.cds_id = c.id
    JOIN hsp h ON h.cds_id = c.id
    WHERE c.gbk_id = ?
      AND h.accession != ''
      AND h.bit_score >= 20
    ORDER BY c.nt_start ASC, h.bit_score DESC
"""

FAMILY_QUERY = """
    SELECT brf.family_id, f.cutoff
    FROM bgc_record_family brf
    JOIN family f ON brf.family_id = f.id
    WHERE brf.record_id = ?
    ORDER BY f.cutoff
"""


# ─── Ordered domain sequence extraction ──────────────────────────────────────

def get_ordered_sequence(cur, gbk_id):
    """
    Return the ordered list of best-scoring Pfam accessions for a BGC,
    one per annotated CDS, sorted by genomic position.
    """
    cur.execute(DOMAIN_SEQ_QUERY, (gbk_id,))
    rows = cur.fetchall()

    # Take the best-scoring domain per CDS (rows already ordered by score DESC)
    seen_cds = {}
    ordered = []
    for row in rows:
        cds_id = row['cds_id']
        if cds_id not in seen_cds:
            seen_cds[cds_id] = True
            ordered.append(row['accession'].split('.')[0])  # strip version

    return ordered


# ─── LCS distance ─────────────────────────────────────────────────────────────

def lcs_length(a, b):
    """Dynamic programming LCS on lists of domain accessions."""
    m, n = len(a), len(b)
    # Use two rows to save memory
    prev = [0] * (n + 1)
    curr = [0] * (n + 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if a[i - 1] == b[j - 1]:
                curr[j] = prev[j - 1] + 1
            else:
                curr[j] = max(prev[j], curr[j - 1])
        prev, curr = curr, [0] * (n + 1)
    return prev[n]


def lcs_distance(a, b):
    """
    Normalized LCS distance: 1 - LCS(a,b) / max(|a|, |b|).
    Returns 1.0 if either sequence is empty.
    """
    if not a or not b:
        return 1.0
    lcs = lcs_length(a, b)
    return 1.0 - lcs / max(len(a), len(b))


# ─── Core pipeline ────────────────────────────────────────────────────────────

def load_bgc_sequences(db_path, bgc_type_filter, family_id=None):
    labels, sequences, metadata = [], [], []
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
            raise ValueError(f"No BGCs found matching '{bgc_type_filter}'")

        print(f"  Found {len(bgc_rows)} BGCs")

        for row in bgc_rows:
            seq = get_ordered_sequence(cur, row['gbk_id'])
            if not seq:
                skipped += 1
                continue

            cur.execute(FAMILY_QUERY, (row['bgc_id'],))
            families = [{'family_id': r['family_id'], 'cutoff': r['cutoff']}
                        for r in cur.fetchall()]

            label = label_from_path(row['gbk_path'])
            labels.append(label)
            sequences.append(seq)
            metadata.append({
                'label':    label,
                'product':  row['product'],
                'organism': row['organism'] or os.path.basename(os.path.dirname(row['gbk_path'])),
                'gbk_path': row['gbk_path'],
                'n_genes':  len(seq),
                'families': families,
                'sequence': seq,                          # accessions
                'sequence_named': [domain_name(d) for d in seq],  # readable names
            })

    if skipped:
        print(f"  Skipped {skipped} BGCs with no annotated domains")

    labels = make_labels_unique(labels)
    for i, m in enumerate(metadata):
        m['label'] = labels[i]

    return labels, sequences, metadata


def build_distance_matrix(labels, sequences):
    n = len(labels)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = lcs_distance(sequences[i], sequences[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist


# ─── Output ───────────────────────────────────────────────────────────────────

def save_outputs(outdir, bgc_type, labels, sequences, metadata, dist_matrix, tree):
    os.makedirs(outdir, exist_ok=True)

    # Newick tree
    newick_path = os.path.join(outdir, f'{bgc_type}_synteny_tree.nwk')
    Phylo.write(tree, newick_path, 'newick')
    print(f"  Tree:             {newick_path}")

    # Ordered domain sequences (human-readable)
    seq_path = os.path.join(outdir, f'{bgc_type}_domain_sequences.tsv')
    with open(seq_path, 'w') as f:
        f.write('BGC\tn_genes\tdomain_sequence\tdomain_accessions\n')
        for m in metadata:
            named   = ' - '.join(m['sequence_named'])
            accs    = ' - '.join(m['sequence'])
            f.write(f"{m['label']}\t{m['n_genes']}\t{named}\t{accs}\n")
    print(f"  Domain sequences: {seq_path}")

    # Distance matrix
    dist_path = os.path.join(outdir, f'{bgc_type}_lcs_distances.tsv')
    with open(dist_path, 'w') as f:
        f.write('\t' + '\t'.join(labels) + '\n')
        for i, label in enumerate(labels):
            row = '\t'.join(f'{dist_matrix[i][j]:.4f}' for j in range(len(labels)))
            f.write(f'{label}\t{row}\n')
    print(f"  Distance matrix:  {dist_path}")

    # Metadata JSON
    meta_path = os.path.join(outdir, f'{bgc_type}_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"  Metadata:         {meta_path}")

    # ASCII tree preview
    print('\n--- Tree preview (ASCII) ---')
    buf = StringIO()
    Phylo.draw_ascii(tree, file=buf)
    lines = buf.getvalue().splitlines()
    for line in lines[:40]:
        print(line)
    if len(lines) > 40:
        print(f'  ... ({len(lines) - 40} more lines — open the .nwk file for the full tree)')


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build a gene-order (synteny) NJ tree of BGCs using LCS distance',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--db',        required=True,
                        help='Path to BiG-SCAPE SQLite database')
    parser.add_argument('--bgc_type',  default='phosphonate',
                        help='BGC product type filter (default: phosphonate)')
    parser.add_argument('--outdir',    required=True,
                        help='Output directory')
    parser.add_argument('--family_id', type=int, default=None,
                        help='Restrict to a single GCF family (e.g. --family_id 2)')
    args = parser.parse_args()

    print(f'Loading BGC domain sequences from: {args.db}')
    labels, sequences, metadata = load_bgc_sequences(
        args.db, args.bgc_type, family_id=args.family_id)

    if len(labels) < 3:
        raise ValueError(f'Need at least 3 BGCs, found {len(labels)}')

    lens = [len(s) for s in sequences]
    print(f'  Gene sequence length: min={min(lens)}, max={max(lens)}, '
          f'mean={sum(lens)/len(lens):.1f}')

    print('Computing pairwise LCS distances...')
    dist_matrix = build_distance_matrix(labels, sequences)

    print('Constructing NJ tree...')
    tree = build_nj_tree(labels, dist_matrix)

    print('Saving outputs...')
    save_outputs(args.outdir, args.bgc_type, labels, sequences,
                 metadata, dist_matrix, tree)

    print('\nDone.')


if __name__ == '__main__':
    main()
