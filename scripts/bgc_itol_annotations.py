#!/usr/bin/env python3
"""
Generate iTOL annotation files for a BGC Pfam NJ tree.

Reads outputs from bgc_pfam_tree.py and produces three annotation files
ready to upload alongside the .nwk file in iTOL (https://itol.embl.de).

Outputs:
  {bgc_type}_itol_gcf.txt       - GCF family membership color strip
  {bgc_type}_itol_domains.txt   - Key pathway domain presence/absence (binary)
  {bgc_type}_itol_domaincount.txt - Total Pfam domain count per BGC (bar chart)

Usage:
    python scripts/bgc_itol_annotations.py \
        --treedir results/bgc_trees/Pantoea \
        --bgc_type phosphonate
"""

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.colors import family_color


# ─── Key pathway domain definitions ──────────────────────────────────────────
# Pfam accessions verified against antiSMASH clusterhmmer output.
# Edit labels/colors here to customise the binary annotation strip.
# Shape codes: 1=circle, 2=triangle, 3=square, 4=star, 5=diamond
KEY_DOMAINS = [
    # (pfam_accession, display_label, hex_color, shape)
    ('PF13714', 'PEP mutase (pepM/aepX)',         '#d62728', 1),  # hallmark — universal
    ('PF00296', 'HEPD / luciferase-like',          '#ff7f0e', 2),  # 2-HEP dioxygenase
    ('PF02775', 'ThDP-binding (Ppd)',               '#2ca02c', 1),  # phosphonopyruvate decarboxylase
    ('PF00266', 'Aminotransferase class V (AepC)',  '#1f77b4', 1),  # 2-aminoethylphosphonate transaminase
    ('PF00682', 'FrbC-like (PmmS)',                '#8c564b', 3),  # phosphonomethylmalate synthase
    ('PF13649', 'Radical SAM',                     '#9467bd', 2),  # radical C–P chemistry
    ('PF08241', 'Methyltransferase',               '#e377c2', 1),  # phosphonate methylation
]


# ─── iTOL file writers ────────────────────────────────────────────────────────

def write_gcf_colorstrip(metadata, outpath, bgc_type):
    """
    DATASET_COLORSTRIP: one color per leaf based on GCF family membership.
    Top 20 families get vivid, evenly-spaced hues. Singletons = gray.
    """
    # Count members per family
    family_counts = {}
    for bgc in metadata:
        fid = bgc['families'][0]['family_id'] if bgc['families'] else None
        if fid is not None:
            family_counts[fid] = family_counts.get(fid, 0) + 1

    # Rank families by size
    ranked = sorted(family_counts, key=lambda x: -family_counts[x])
    family_rank = {fid: i for i, fid in enumerate(ranked)}
    n_top = min(20, len(ranked))

    with open(outpath, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tGCF Family ({bgc_type})\n')
        f.write('COLOR\t#777777\n')
        f.write('LEGEND_TITLE\tGCF Family\n')

        # Legend for top families
        legend_shapes = '\t'.join(['1'] * n_top + ['1'])
        legend_colors = '\t'.join(
            [family_color(fid, family_rank[fid], n_top) for fid in ranked[:n_top]]
            + ['#cccccc']
        )
        legend_labels = '\t'.join(
            [f'GCF-{fid} ({family_counts[fid]} BGCs)' for fid in ranked[:n_top]]
            + ['Singleton / unclustered']
        )
        f.write(f'LEGEND_SHAPES\t{legend_shapes}\n')
        f.write(f'LEGEND_COLORS\t{legend_colors}\n')
        f.write(f'LEGEND_LABELS\t{legend_labels}\n')
        f.write('DATA\n')

        for bgc in metadata:
            fid = bgc['families'][0]['family_id'] if bgc['families'] else None
            rank = family_rank.get(fid, len(ranked))
            color = family_color(fid, rank, n_top)
            label = f"GCF-{fid}" if fid is not None else "Singleton"
            f.write(f"{bgc['label']}\t{color}\t{label}\n")

    print(f"  GCF color strip:  {outpath}  ({len(set(family_counts))} families, "
          f"{sum(1 for b in metadata if not b['families'])} singletons)")


def write_domain_binary(metadata, domain_matrix_path, outpath, bgc_type):
    """
    DATASET_BINARY: presence/absence squares for each key pathway domain.
    Reads from domain_matrix.tsv (Pfam tree) or derives from metadata sequences
    (synteny tree) if the matrix file is absent.
    """
    # Determine which key domains are present across all BGCs
    if domain_matrix_path and os.path.exists(domain_matrix_path):
        # Read domain matrix header to get column indices
        with open(domain_matrix_path) as f:
            header = f.readline().rstrip('\n').split('\t')[1:]
        domain_index = {d: i for i, d in enumerate(header)}
        active = [(acc, lbl, col, shp) for acc, lbl, col, shp in KEY_DOMAINS
                  if acc in domain_index]
        missing = [acc for acc, *_ in KEY_DOMAINS if acc not in domain_index]
        if missing:
            print(f"  Note: domains not found in matrix (will be skipped): {missing}")
        bgc_presence = {}
        label_set = {b['label'] for b in metadata}
        with open(domain_matrix_path) as f:
            f.readline()
            for line in f:
                parts = line.rstrip('\n').split('\t')
                lbl = parts[0]
                if lbl in label_set:
                    vals = parts[1:]
                    bgc_presence[lbl] = {acc: int(vals[domain_index[acc]]) for acc, *_ in active}
    else:
        # Derive presence/absence from sequences stored in metadata
        all_domains = set()
        for bgc in metadata:
            all_domains.update(bgc.get('sequence', []))
        active = [(acc, lbl, col, shp) for acc, lbl, col, shp in KEY_DOMAINS
                  if acc in all_domains]
        missing = [acc for acc, *_ in KEY_DOMAINS if acc not in all_domains]
        if missing:
            print(f"  Note: domains not found in sequences (will be skipped): {missing}")
        bgc_presence = {
            bgc['label']: {acc: (1 if acc in bgc.get('sequence', []) else 0)
                           for acc, *_ in active}
            for bgc in metadata
        }

    with open(outpath, 'w') as f:
        f.write('DATASET_BINARY\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\tKey domains ({bgc_type})\n')
        f.write('COLOR\t#333333\n')

        f.write('FIELD_SHAPES\t' + '\t'.join(str(shp) for _, _, _, shp in active) + '\n')
        f.write('FIELD_LABELS\t' + '\t'.join(lbl for _, lbl, _, _ in active) + '\n')
        f.write('FIELD_COLORS\t' + '\t'.join(col for _, _, col, _ in active) + '\n')

        f.write('LEGEND_TITLE\tPathway domains\n')
        f.write('LEGEND_SHAPES\t' + '\t'.join(str(shp) for _, _, _, shp in active) + '\n')
        f.write('LEGEND_COLORS\t' + '\t'.join(col for _, _, col, _ in active) + '\n')
        f.write('LEGEND_LABELS\t' + '\t'.join(lbl for _, lbl, _, _ in active) + '\n')

        f.write('DATA\n')
        for bgc in metadata:
            presence = bgc_presence.get(bgc['label'], {acc: 0 for acc, *_ in active})
            vals = '\t'.join(str(presence.get(acc, 0)) for acc, *_ in active)
            f.write(f"{bgc['label']}\t{vals}\n")

    print(f"  Domain binary:    {outpath}  ({len(active)} domains shown)")


def write_domain_count_bar(metadata, outpath, bgc_type, count_field='n_domains'):
    """
    DATASET_SIMPLEBAR: total domain/gene count per BGC.
    """
    label_str = 'Gene count' if count_field == 'n_genes' else 'Pfam domain count'
    with open(outpath, 'w') as f:
        f.write('DATASET_SIMPLEBAR\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{label_str} ({bgc_type})\n')
        f.write('COLOR\t#5b5ea6\n')
        f.write('WIDTH\t200\n')
        f.write('SHOW_INTERNAL\t0\n')
        f.write('DATA\n')
        for bgc in metadata:
            f.write(f"{bgc['label']}\t{bgc[count_field]}\n")

    counts = [b[count_field] for b in metadata]
    print(f"  Domain count bar: {outpath}  "
          f"(range {min(counts)}–{max(counts)}, mean {sum(counts)/len(counts):.1f})")


# ─── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate iTOL annotation files for a BGC Pfam NJ tree',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--treedir',  required=True,
                        help='Directory containing bgc_pfam_tree.py outputs')
    parser.add_argument('--bgc_type', default='phosphonate',
                        help='BGC type prefix used in filenames (default: phosphonate)')
    args = parser.parse_args()

    d = args.treedir
    t = args.bgc_type

    meta_path   = os.path.join(d, f'{t}_metadata.json')
    matrix_path = os.path.join(d, f'{t}_domain_matrix.tsv')  # optional

    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Required input not found: {meta_path}\n"
                                "Run bgc_pfam_tree.py or bgc_synteny_tree.py first.")

    with open(meta_path) as f:
        metadata = json.load(f)

    print(f"Generating iTOL annotations for {len(metadata)} BGCs...")

    write_gcf_colorstrip(
        metadata,
        os.path.join(d, f'{t}_itol_gcf.txt'),
        t
    )
    write_domain_binary(
        metadata,
        matrix_path if os.path.exists(matrix_path) else None,
        os.path.join(d, f'{t}_itol_domains.txt'),
        t
    )
    # Use n_genes for synteny tree, n_domains for Pfam tree
    count_field = 'n_genes' if 'n_genes' in metadata[0] else 'n_domains'
    write_domain_count_bar(
        metadata,
        os.path.join(d, f'{t}_itol_domaincount.txt'),
        t,
        count_field=count_field
    )

    print(f"\nDone. Upload these files to iTOL alongside the .nwk tree:")
    print(f"  {t}_pfam_tree.nwk")
    print(f"  {t}_itol_gcf.txt")
    print(f"  {t}_itol_domains.txt")
    print(f"  {t}_itol_domaincount.txt")


if __name__ == '__main__':
    main()
