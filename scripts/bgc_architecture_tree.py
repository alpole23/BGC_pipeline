#!/usr/bin/env python3
"""
Build a tree of unique BGC domain architectures within a GCF.

BGCs are grouped by their exact domain composition (unordered multiset).
Identical BGCs collapse to a single leaf annotated with a count. Pairwise
distances use generalized Jaccard on domain multisets, which handles both
presence/absence and copy number variation.

    distance(A, B) = 1 - sum(min(A[d], B[d])) / sum(max(A[d], B[d]))

This approach:
  - Is orientation-independent (no strand issues)
  - Collapses redundant BGCs so each leaf = one unique architecture
  - Preserves copy number (2x ATPgrasp ≠ 1x ATPgrasp)
  - Separates common core architectures from rare/novel variants

Usage:
    python scripts/bgc_architecture_tree.py \\
        --db results/bigscape_results/Pantoea/Pantoea.db \\
        --bgc_type phosphonate \\
        --family_id 2 \\
        --outdir results/bgc_trees/Pantoea/GCF2_arch
"""

import argparse
import colorsys
import hashlib
import sqlite3
import os
import sys
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import Phylo
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import DOMAIN_NAMES
from utils.tree_building import build_nj_tree
from utils.colors import genus_color as _genus_color

sys.setrecursionlimit(10000)


def domain_name(acc):
    base = acc.split('.')[0]
    return DOMAIN_NAMES.get(base, base)


# ─── Database queries ─────────────────────────────────────────────────────────

BGC_QUERY = """
    SELECT br.id AS bgc_id, br.product,
           g.id AS gbk_id, g.path AS gbk_path, g.organism
    FROM bgc_record br
    JOIN gbk g ON br.gbk_id = g.id
    WHERE LOWER(br.product) LIKE ?
    ORDER BY g.path
"""

BGC_QUERY_FAMILY = """
    SELECT br.id AS bgc_id, br.product,
           g.id AS gbk_id, g.path AS gbk_path, g.organism
    FROM bgc_record br
    JOIN gbk g ON br.gbk_id = g.id
    JOIN bgc_record_family brf ON br.id = brf.record_id
    WHERE LOWER(br.product) LIKE ?
      AND brf.family_id = ?
    ORDER BY g.path
"""

DOMAIN_QUERY = """
    SELECT h.accession, h.bit_score
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


# ─── Domain extraction ────────────────────────────────────────────────────────

def get_domain_counter(cur, gbk_id):
    """
    Return a Counter of domain accessions for a BGC.
    Takes the best-scoring domain per CDS; counts multiple copies separately.
    """
    cur.execute("""
        SELECT c.id as cds_id, c.nt_start, h.accession, h.bit_score
        FROM cds c
        JOIN scanned_cds sc ON sc.cds_id = c.id
        JOIN hsp h ON h.cds_id = c.id
        WHERE c.gbk_id = ?
          AND h.accession != ''
          AND h.bit_score >= 20
        ORDER BY c.nt_start ASC, h.bit_score DESC
    """, (gbk_id,))
    rows = cur.fetchall()

    seen_cds = {}
    counter = Counter()
    for row in rows:
        cds_id = row[0]
        if cds_id not in seen_cds:
            seen_cds[cds_id] = True
            acc = row[2].split('.')[0]
            counter[acc] += 1

    return counter


# ─── Distance ─────────────────────────────────────────────────────────────────

def multiset_jaccard_distance(counter_a, counter_b):
    """
    Generalized Jaccard distance for domain multisets.
    Handles copy number: 2x ATPgrasp ≠ 1x ATPgrasp.
    Returns 0.0 for identical architectures, 1.0 for completely disjoint.
    """
    all_keys = set(counter_a) | set(counter_b)
    intersection = sum(min(counter_a.get(k, 0), counter_b.get(k, 0)) for k in all_keys)
    union        = sum(max(counter_a.get(k, 0), counter_b.get(k, 0)) for k in all_keys)
    return 1.0 - intersection / union if union > 0 else 1.0


# ─── Core pipeline ────────────────────────────────────────────────────────────

def load_and_group(db_path, bgc_type_filter, family_id=None):
    """
    Load all BGCs, extract domain counters, and group by identical architecture.

    Returns:
        architectures : list of dicts, one per unique architecture, sorted by count desc
    """
    skipped = 0
    groups = defaultdict(lambda: {'counter': None, 'genomes': [], 'organisms': []})

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

        print(f"  Found {len(bgc_rows)} BGCs — grouping by domain architecture...")

        for row in bgc_rows:
            counter = get_domain_counter(cur, row['gbk_id'])
            if not counter:
                skipped += 1
                continue

            key = tuple(sorted(counter.items()))  # canonical hashable key
            if groups[key]['counter'] is None:
                groups[key]['counter'] = counter
            genome_name = os.path.basename(os.path.dirname(row['gbk_path']))
            groups[key]['genomes'].append(genome_name)
            groups[key]['organisms'].append(row['organism'] or genome_name)

    if skipped:
        print(f"  Skipped {skipped} BGCs with no domain annotations")

    # Sort by count descending
    sorted_groups = sorted(groups.values(), key=lambda g: -len(g['genomes']))

    # Build architecture records with readable labels
    architectures = []
    for i, grp in enumerate(sorted_groups, start=1):
        count    = len(grp['genomes'])
        counter  = grp['counter']
        label    = f"arch_{i:03d}_n{count}"

        # Dominant organism (most frequent genus)
        genus_counts = Counter(
            o.split()[0] if o else 'Unknown' for o in grp['organisms']
        )
        dominant_genus = genus_counts.most_common(1)[0][0]

        # Human-readable domain summary (sorted by count desc, then name)
        domain_parts = []
        for acc, cnt in sorted(counter.items(), key=lambda x: (-x[1], x[0])):
            name = domain_name(acc)
            domain_parts.append(f"{name}x{cnt}" if cnt > 1 else name)
        domain_str = ' | '.join(domain_parts)

        architectures.append({
            'label':          label,
            'rank':           i,
            'count':          count,
            'counter':        dict(counter),
            'domain_str':     domain_str,
            'n_domain_types': len(counter),
            'total_domains':  sum(counter.values()),
            'dominant_genus': dominant_genus,
            'genomes':        grp['genomes'],
            'organisms':      grp['organisms'],
        })

    print(f"  {len(bgc_rows) - skipped} BGCs → {len(architectures)} unique architectures")
    return architectures


def build_distance_matrix(architectures):
    n = len(architectures)
    counters = [Counter(a['counter']) for a in architectures]
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = multiset_jaccard_distance(counters[i], counters[j])
            dist[i, j] = d
            dist[j, i] = d
    return dist


# ─── Output ───────────────────────────────────────────────────────────────────

def save_outputs(outdir, bgc_type, architectures, dist_matrix, tree):
    os.makedirs(outdir, exist_ok=True)
    labels = [a['label'] for a in architectures]

    # Newick tree
    newick_path = os.path.join(outdir, f'{bgc_type}_arch_tree.nwk')
    Phylo.write(tree, newick_path, 'newick')
    print(f"  Tree:              {newick_path}")

    # Architecture summary TSV
    arch_path = os.path.join(outdir, f'{bgc_type}_architectures.tsv')
    with open(arch_path, 'w') as f:
        f.write('label\tcount\tn_domain_types\ttotal_domains\t'
                'dominant_genus\tdomain_composition\texample_genomes\n')
        for a in architectures:
            examples = ', '.join(a['genomes'][:3])
            if len(a['genomes']) > 3:
                examples += f' (+{len(a["genomes"]) - 3} more)'
            f.write(f"{a['label']}\t{a['count']}\t{a['n_domain_types']}\t"
                    f"{a['total_domains']}\t{a['dominant_genus']}\t"
                    f"{a['domain_str']}\t{examples}\n")
    print(f"  Architecture TSV:  {arch_path}")

    # Distance matrix
    dist_path = os.path.join(outdir, f'{bgc_type}_arch_distances.tsv')
    with open(dist_path, 'w') as f:
        f.write('\t' + '\t'.join(labels) + '\n')
        for i, label in enumerate(labels):
            row = '\t'.join(f'{dist_matrix[i][j]:.4f}' for j in range(len(labels)))
            f.write(f'{label}\t{row}\n')
    print(f"  Distance matrix:   {dist_path}")

    # Metadata JSON
    meta_path = os.path.join(outdir, f'{bgc_type}_metadata.json')
    with open(meta_path, 'w') as f:
        json.dump([{k: v for k, v in a.items() if k != 'counter'}
                   for a in architectures], f, indent=2)
    print(f"  Metadata:          {meta_path}")

    # iTOL files
    _write_itol_count_bar(architectures, outdir, bgc_type)
    _write_itol_domain_binary(architectures, outdir, bgc_type)
    _write_itol_genus_colorstrip(architectures, outdir, bgc_type)
    _write_itol_arch_label(architectures, outdir, bgc_type)
    _write_itol_genome_list(architectures, outdir, bgc_type)

    # ASCII preview
    print('\n--- Tree preview (ASCII) ---')
    buf = StringIO()
    Phylo.draw_ascii(tree, file=buf)
    lines = buf.getvalue().splitlines()
    for line in lines[:50]:
        print(line)
    if len(lines) > 50:
        print(f'  ... ({len(lines) - 50} more lines)')


# ─── iTOL writers ─────────────────────────────────────────────────────────────

def _write_itol_count_bar(architectures, outdir, bgc_type):
    path = os.path.join(outdir, f'{bgc_type}_itol_count.txt')
    with open(path, 'w') as f:
        f.write('DATASET_SIMPLEBAR\nSEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tBGC count\n')
        f.write('COLOR\t#2c7bb6\nWIDTH\t300\nSHOW_INTERNAL\t0\nDATA\n')
        for a in architectures:
            f.write(f"{a['label']}\t{a['count']}\n")
    print(f"  iTOL count bar:    {path}")


KEY_DOMAINS_ITOL = [
    ('PF13714', 'PEP_mutase',    '#d62728', 1),
    ('PF00296', 'HEPD',          '#ff7f0e', 2),
    ('PF02775', 'ThDP_C (Ppd)',  '#2ca02c', 1),
    ('PF00266', 'Aminotrans_V',  '#1f77b4', 1),
    ('PF00682', 'HMGL',          '#8c564b', 3),
    ('PF13649', 'Radical_SAM',   '#9467bd', 2),
    ('PF08241', 'Methyltransf',  '#e377c2', 1),
    ('PF00330', 'Aconitase',     '#bcbd22', 1),
    ('PF00694', 'Aconitase_C',   '#17becf', 1),
]


def _write_itol_domain_binary(architectures, outdir, bgc_type):
    all_domains = set()
    for a in architectures:
        all_domains.update(a['counter'].keys())
    active = [(acc, lbl, col, shp) for acc, lbl, col, shp in KEY_DOMAINS_ITOL
              if acc in all_domains]

    path = os.path.join(outdir, f'{bgc_type}_itol_domains.txt')
    with open(path, 'w') as f:
        f.write('DATASET_BINARY\nSEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tKey domains\nCOLOR\t#333333\n')
        f.write('FIELD_SHAPES\t'  + '\t'.join(str(s) for _, _, _, s in active) + '\n')
        f.write('FIELD_LABELS\t'  + '\t'.join(l for _, l, _, _ in active) + '\n')
        f.write('FIELD_COLORS\t'  + '\t'.join(c for _, _, c, _ in active) + '\n')
        f.write('LEGEND_TITLE\tPathway domains\n')
        f.write('LEGEND_SHAPES\t' + '\t'.join(str(s) for _, _, _, s in active) + '\n')
        f.write('LEGEND_COLORS\t' + '\t'.join(c for _, _, c, _ in active) + '\n')
        f.write('LEGEND_LABELS\t' + '\t'.join(l for _, l, _, _ in active) + '\n')
        f.write('DATA\n')
        for a in architectures:
            vals = '\t'.join('1' if acc in a['counter'] else '0' for acc, *_ in active)
            f.write(f"{a['label']}\t{vals}\n")
    print(f"  iTOL domain binary:{path}")


def _write_itol_genus_colorstrip(architectures, outdir, bgc_type):
    genus_counts = Counter(a['dominant_genus'] for a in architectures)
    ranked = [g for g, _ in genus_counts.most_common()]
    genus_rank = {g: i for i, g in enumerate(ranked)}
    n_top = min(12, len(ranked))

    path = os.path.join(outdir, f'{bgc_type}_itol_genus.txt')
    with open(path, 'w') as f:
        f.write('DATASET_COLORSTRIP\nSEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tDominant genus\nCOLOR\t#777777\n')
        f.write('LEGEND_TITLE\tGenus\n')
        f.write('LEGEND_SHAPES\t' + '\t'.join(['1'] * n_top) + '\n')
        f.write('LEGEND_COLORS\t' + '\t'.join(
            _genus_color(g, genus_rank[g], n_top) for g in ranked[:n_top]) + '\n')
        f.write('LEGEND_LABELS\t' + '\t'.join(ranked[:n_top]) + '\n')
        f.write('DATA\n')
        for a in architectures:
            g = a['dominant_genus']
            color = _genus_color(g, genus_rank[g], n_top)
            f.write(f"{a['label']}\t{color}\t{g}\n")
    print(f"  iTOL genus strip:  {path}")


def _write_itol_arch_label(architectures, outdir, bgc_type):
    """DATASET_TEXT: show the domain composition string next to each leaf."""
    path = os.path.join(outdir, f'{bgc_type}_itol_archlabel.txt')
    with open(path, 'w') as f:
        f.write('DATASET_TEXT\nSEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tDomain architecture\nCOLOR\t#333333\n')
        f.write('SHOW_INTERNAL\t0\nSIZE_FACTOR\t0.8\nDATA\n')
        for a in architectures:
            label_str = a['domain_str']
            if len(label_str) > 120:
                label_str = label_str[:117] + '...'
            f.write(f"{a['label']}\t{label_str}\t-1\t#333333\tnormal\t1\n")
    print(f"  iTOL arch labels:  {path}")


def _write_itol_genome_list(architectures, outdir, bgc_type):
    """
    DATASET_TEXT: genome names for each architecture leaf.
    - n=1      : show full genome name
    - n=2–10   : show all genome names pipe-separated
    - n>10     : show first 5 names + '(+N more)'
    """
    path = os.path.join(outdir, f'{bgc_type}_itol_genomes.txt')
    with open(path, 'w') as f:
        f.write('DATASET_TEXT\nSEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tGenomes\nCOLOR\t#333333\n')
        f.write('SHOW_INTERNAL\t0\nSIZE_FACTOR\t0.8\nDATA\n')
        for a in architectures:
            genomes = a['genomes']
            n = len(genomes)
            if n <= 10:
                text = ' | '.join(genomes)
            else:
                text = ' | '.join(genomes[:5]) + f' (+{n - 5} more)'
            f.write(f"{a['label']}\t{text}\t-1\t#555555\tnormal\t1\n")
    print(f"  iTOL genome list:  {path}")


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build a tree of unique BGC architectures within a GCF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--db',        required=True)
    parser.add_argument('--bgc_type',  default='phosphonate')
    parser.add_argument('--outdir',    required=True)
    parser.add_argument('--family_id', type=int, default=None,
                        help='Restrict to one GCF (e.g. --family_id 2)')
    args = parser.parse_args()

    print(f'Loading BGCs from: {args.db}')
    architectures = load_and_group(args.db, args.bgc_type, family_id=args.family_id)

    if len(architectures) < 3:
        raise ValueError(f'Need at least 3 unique architectures, '
                         f'found {len(architectures)}')

    print(f'Computing pairwise distances ({len(architectures)} architectures)...')
    dist_matrix = build_distance_matrix(architectures)

    print('Constructing NJ tree...')
    labels = [a['label'] for a in architectures]
    tree = build_nj_tree(labels, dist_matrix)

    print('Saving outputs...')
    save_outputs(args.outdir, args.bgc_type, architectures, dist_matrix, tree)

    print(f'\nDone. Upload to iTOL:')
    print(f'  {args.bgc_type}_arch_tree.nwk')
    print(f'  {args.bgc_type}_itol_count.txt')
    print(f'  {args.bgc_type}_itol_domains.txt')
    print(f'  {args.bgc_type}_itol_genus.txt')
    print(f'  {args.bgc_type}_itol_archlabel.txt')


if __name__ == '__main__':
    main()
