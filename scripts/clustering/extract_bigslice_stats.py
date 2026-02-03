#!/usr/bin/env python3
"""Extract statistics from BiG-SLiCE clustering results."""

import argparse
import json
import os
import sqlite3


def extract_stats(bigslice_dir, output_file):
    # Path to BiG-SLiCE database
    db_path = os.path.join(bigslice_dir, "bigslice_output", "result", "data.db")

    stats = {}

    if not os.path.exists(db_path):
        stats["error"] = f"Database not found: {db_path}"
        with open(output_file, "w") as f:
            json.dump(stats, f, indent=2)
        return

    with sqlite3.connect(db_path) as con:
        cur = con.cursor()

        # Get run information
        run_info = cur.execute("SELECT id, threshold FROM clustering WHERE run_id=1").fetchone()
        if run_info:
            clustering_id, threshold = run_info
            stats["clustering_id"] = clustering_id
            stats["threshold"] = float(threshold)
        else:
            # No clustering data
            stats["error"] = "No clustering data found"
            with open(output_file, "w") as f:
                json.dump(stats, f, indent=2)
            return

        # Basic clustering statistics
        result = cur.execute(
            """SELECT
                COUNT(DISTINCT bgc.id) as total_bgcs,
                COUNT(DISTINCT gcf.id) as total_gcfs,
                AVG(membership_value) as avg_distance
            FROM bgc, gcf, gcf_membership
            WHERE gcf_membership.bgc_id=bgc.id
              AND gcf_membership.gcf_id=gcf.id
              AND gcf.clustering_id=?""",
            (clustering_id,)).fetchone()

        stats["total_bgcs"] = result[0] if result[0] else 0
        stats["total_gcfs"] = result[1] if result[1] else 0
        stats["avg_distance_to_gcf"] = round(result[2], 2) if result[2] else 0

        # Threshold and assignment statistics
        result = cur.execute(
            """SELECT
                COUNT(DISTINCT bgc_id),
                COUNT(DISTINCT gcf_id)
            FROM gcf_membership
            WHERE gcf_id IN (SELECT id FROM gcf WHERE clustering_id=?)
              AND rank=0
              AND membership_value <= ?""",
            (clustering_id, threshold)).fetchone()

        stats["bgcs_assigned"] = result[0] if result[0] else 0
        stats["gcfs_with_members"] = result[1] if result[1] else 0

        stats["bgcs_not_assigned"] = stats["total_bgcs"] - stats["bgcs_assigned"]

        # BGCs per GCF distribution
        gcf_distribution = cur.execute(
            """SELECT
                gcf.id as gcf_id,
                COUNT(bgc.id) as bgc_count
            FROM bgc, gcf, gcf_membership
            WHERE gcf_membership.bgc_id=bgc.id
              AND gcf_membership.gcf_id=gcf.id
              AND gcf.clustering_id=?
              AND gcf_membership.rank=0
            GROUP BY gcf.id
            ORDER BY bgc_count DESC""",
            (clustering_id,)).fetchall()

        stats["gcf_distribution"] = [
            {"gcf_id": row[0], "bgc_count": row[1]}
            for row in gcf_distribution
        ]

        # Calculate GCF distribution summary
        if gcf_distribution:
            bgc_counts = [row[1] for row in gcf_distribution]
            stats["min_bgcs_per_gcf"] = min(bgc_counts)
            stats["max_bgcs_per_gcf"] = max(bgc_counts)
            stats["avg_bgcs_per_gcf"] = round(sum(bgc_counts) / len(bgc_counts), 2)
            stats["singleton_gcfs"] = sum(1 for count in bgc_counts if count == 1)
        else:
            stats["min_bgcs_per_gcf"] = 0
            stats["max_bgcs_per_gcf"] = 0
            stats["avg_bgcs_per_gcf"] = 0
            stats["singleton_gcfs"] = 0

        # Fragmentation statistics
        frag_stats = cur.execute(
            """SELECT
                bgc.on_contig_edge,
                COUNT(*) as count,
                AVG(membership_value) as avg_distance
            FROM bgc, gcf, gcf_membership
            WHERE gcf_membership.bgc_id=bgc.id
              AND gcf_membership.gcf_id=gcf.id
              AND gcf.clustering_id=?
              AND gcf_membership.rank=0
            GROUP BY bgc.on_contig_edge""",
            (clustering_id,)).fetchall()

        stats["fragmentation"] = {}
        for row in frag_stats:
            edge_status = row[0]
            if edge_status == 0:
                frag_type = "complete"
            elif edge_status == 1:
                frag_type = "fragmented"
            else:
                frag_type = "unknown"

            stats["fragmentation"][frag_type] = {
                "count": row[1],
                "avg_distance": round(row[2], 2) if row[2] else 0
            }

    # Write statistics to JSON file
    with open(output_file, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"BiG-SLiCE statistics extracted successfully")
    print(f"Total BGCs: {stats['total_bgcs']}")
    print(f"Total GCFs: {stats['total_gcfs']}")
    print(f"BGCs assigned (threshold={stats['threshold']}): {stats['bgcs_assigned']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract BiG-SLiCE clustering statistics")
    parser.add_argument("bigslice_dir", help="Path to BiG-SLiCE output directory")
    parser.add_argument("output_file", help="Path to output JSON file")
    args = parser.parse_args()

    extract_stats(args.bigslice_dir, args.output_file)
