#!/usr/bin/env python3
"""
Generate rarefaction curves from BiG-SCAPE clustering results.

Shows how the number of unique Gene Cluster Families (GCFs) discovered
increases as more genomes are sampled.
"""

import argparse
import sqlite3
import os
import re
from pathlib import Path
from collections import defaultdict
import random

import numpy as np
import matplotlib.pyplot as plt


def extract_genome_name(path):
    """Extract genome name from BiG-SCAPE GBK path."""
    # Path format: .../antismash_input/Genome_Name/contig.regionXXX.gbk
    parts = path.split('/')
    for i, part in enumerate(parts):
        if part == 'antismash_input' and i + 1 < len(parts):
            return parts[i + 1]
    # Fallback: use parent directory name
    return Path(path).parent.name


def load_data_from_db(db_path, bgc_type=None):
    """
    Load genome -> GCF mapping from BiG-SCAPE database.

    Returns:
        dict: {genome_name: set(gcf_ids)}
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Query to get genome path, BGC product, and GCF assignment
    if bgc_type:
        query = """
            SELECT g.path, b.product, bf.family_id
            FROM gbk g
            JOIN bgc_record b ON g.id = b.gbk_id
            LEFT JOIN bgc_record_family bf ON b.id = bf.record_id
            WHERE b.product LIKE ?
        """
        cursor.execute(query, (f'%{bgc_type}%',))
    else:
        query = """
            SELECT g.path, b.product, bf.family_id
            FROM gbk g
            JOIN bgc_record b ON g.id = b.gbk_id
            LEFT JOIN bgc_record_family bf ON b.id = bf.record_id
        """
        cursor.execute(query)

    genome_gcfs = defaultdict(set)
    bgc_counts = defaultdict(int)

    for path, product, family_id in cursor.fetchall():
        genome = extract_genome_name(path)
        if family_id is not None:
            genome_gcfs[genome].add(family_id)
        bgc_counts[genome] += 1

    conn.close()
    return dict(genome_gcfs), dict(bgc_counts)


def calculate_rarefaction(genome_gcfs, n_iterations=100, n_points=50):
    """
    Calculate rarefaction curve with bootstrap confidence intervals.

    Args:
        genome_gcfs: dict mapping genome names to sets of GCF IDs
        n_iterations: number of bootstrap iterations
        n_points: number of points on the curve

    Returns:
        x_values: genome counts
        mean_gcfs: mean GCF count at each point
        lower_ci: 2.5th percentile
        upper_ci: 97.5th percentile
    """
    genomes = list(genome_gcfs.keys())
    n_genomes = len(genomes)

    if n_genomes == 0:
        return [], [], [], []

    # Sample points along the curve
    x_values = np.linspace(1, n_genomes, min(n_points, n_genomes)).astype(int)
    x_values = sorted(set(x_values))  # Remove duplicates

    # Bootstrap iterations
    all_curves = []

    for _ in range(n_iterations):
        # Shuffle genome order
        shuffled = genomes.copy()
        random.shuffle(shuffled)

        # Accumulate GCFs
        seen_gcfs = set()
        curve = []

        for i, genome in enumerate(shuffled, 1):
            seen_gcfs.update(genome_gcfs[genome])
            if i in x_values:
                curve.append(len(seen_gcfs))

        all_curves.append(curve)

    # Calculate statistics
    all_curves = np.array(all_curves)
    mean_gcfs = np.mean(all_curves, axis=0)
    lower_ci = np.percentile(all_curves, 2.5, axis=0)
    upper_ci = np.percentile(all_curves, 97.5, axis=0)

    return x_values, mean_gcfs, lower_ci, upper_ci


def get_bgc_types(db_path, min_count=50):
    """Get BGC types with at least min_count occurrences."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT product, COUNT(*) as count
        FROM bgc_record
        GROUP BY product
        HAVING count >= ?
        ORDER BY count DESC
    """, (min_count,))

    types = [(row[0], row[1]) for row in cursor.fetchall()]
    conn.close()
    return types


def plot_rarefaction(results, outdir, taxon):
    """
    Plot rarefaction curves.

    Args:
        results: dict of {label: (x, mean, lower, upper, color)}
        outdir: output directory
        taxon: taxon name for title
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    for label, (x, mean, lower, upper, color) in results.items():
        if len(x) == 0:
            continue
        ax.plot(x, mean, color=color, linewidth=2, label=label)
        ax.fill_between(x, lower, upper, color=color, alpha=0.2)

    ax.set_xlabel('Number of Genomes Sampled', fontsize=12)
    ax.set_ylabel('Unique Gene Cluster Families (GCFs)', fontsize=12)
    ax.set_title(f'GCF Rarefaction Curve - {taxon}', fontsize=14)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add annotation about curve interpretation
    ax.text(0.02, 0.98,
            'Plateau = diversity saturated\nRising = more diversity to discover',
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    output_path = Path(outdir) / 'rarefaction_curve.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Saved rarefaction curve to {output_path}")
    return str(output_path)


def main():
    parser = argparse.ArgumentParser(description='Generate GCF rarefaction curves')
    parser.add_argument('--db', required=True, help='Path to BiG-SCAPE SQLite database')
    parser.add_argument('--outdir', default='.', help='Output directory')
    parser.add_argument('--taxon', default='Unknown', help='Taxon name for title')
    parser.add_argument('--iterations', type=int, default=100, help='Bootstrap iterations')
    parser.add_argument('--types', nargs='*', help='Specific BGC types to plot (default: top 5 + phosphonate)')
    args = parser.parse_args()

    if not os.path.exists(args.db):
        print(f"Error: Database not found: {args.db}")
        return 1

    os.makedirs(args.outdir, exist_ok=True)

    # Colors for different BGC types
    colors = plt.cm.tab10.colors

    results = {}

    # All BGCs combined
    print("Calculating rarefaction for all BGC types...")
    genome_gcfs, bgc_counts = load_data_from_db(args.db)
    n_genomes = len(genome_gcfs)
    n_gcfs = len(set().union(*genome_gcfs.values())) if genome_gcfs else 0
    print(f"  Found {n_genomes} genomes with {n_gcfs} unique GCFs")

    x, mean, lower, upper = calculate_rarefaction(genome_gcfs, args.iterations)
    results['All BGC types'] = (x, mean, lower, upper, 'black')

    # Determine which BGC types to plot
    if args.types:
        types_to_plot = [(t, 0) for t in args.types]
    else:
        # Get top BGC types
        all_types = get_bgc_types(args.db, min_count=50)
        # Always include phosphonate if present
        types_to_plot = all_types[:5]
        phosphonate_types = [t for t in all_types if 'phosphonate' in t[0].lower()]
        if phosphonate_types and phosphonate_types[0] not in types_to_plot:
            types_to_plot.append(phosphonate_types[0])

    # Calculate rarefaction for each BGC type
    for i, (bgc_type, count) in enumerate(types_to_plot):
        print(f"Calculating rarefaction for {bgc_type}...")
        genome_gcfs, _ = load_data_from_db(args.db, bgc_type)
        n_genomes = len(genome_gcfs)
        n_gcfs = len(set().union(*genome_gcfs.values())) if genome_gcfs else 0
        print(f"  Found {n_genomes} genomes with {n_gcfs} unique GCFs")

        if n_genomes > 0:
            x, mean, lower, upper = calculate_rarefaction(genome_gcfs, args.iterations)
            # Simplify label
            label = bgc_type.split('.')[0] if '.' in bgc_type else bgc_type
            results[label] = (x, mean, lower, upper, colors[i % len(colors)])

    # Generate plot
    plot_rarefaction(results, args.outdir, args.taxon)

    # Print summary statistics
    print("\n=== Rarefaction Summary ===")
    for label, (x, mean, lower, upper, _) in results.items():
        if len(x) > 0:
            final_gcfs = mean[-1]
            # Estimate saturation: compare last 10% increase
            if len(mean) >= 10:
                early_rate = (mean[len(mean)//10] - mean[0]) / (x[len(x)//10] - x[0]) if x[len(x)//10] != x[0] else 0
                late_rate = (mean[-1] - mean[-len(mean)//10]) / (x[-1] - x[-len(x)//10]) if x[-1] != x[-len(x)//10] else 0
                saturation = 1 - (late_rate / early_rate) if early_rate > 0 else 1
                print(f"{label}: {final_gcfs:.0f} GCFs at {x[-1]} genomes (saturation: {saturation*100:.1f}%)")
            else:
                print(f"{label}: {final_gcfs:.0f} GCFs at {x[-1]} genomes")

    return 0


if __name__ == '__main__':
    exit(main())
