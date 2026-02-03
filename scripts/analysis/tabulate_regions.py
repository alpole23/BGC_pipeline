#!/usr/bin/env python3
"""Given a bunch of antismash results, tabulate BGC regions."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.constants import KCB_THRESHOLDS


def parse_json(path):
    result_list = []
    with Path.open(path) as f:
        data = json.load(f)
    # Track record index (1-based) for antiSMASH-style region naming
    for record_index, record in enumerate(data["records"], 1):
        if not record["areas"]:
            continue
        regions = [feat for feat in record["features"] if feat["type"] == "region"]
        # Add knownclusterblast if present
        try:
            knownblast = record["modules"]["antismash.modules.clusterblast"][
                "knowncluster"
            ]["results"]
        except (KeyError, TypeError, AttributeError):
            knownblast = None
        for i, region in enumerate(regions):
            # Extract all numbers from location (handles join locations)
            numbers = re.findall(r"\d+", region["location"])
            # For join locations, use min and max; for regular locations, first and last
            start = min(numbers) if len(numbers) > 2 else numbers[0]
            end = max(numbers) if len(numbers) > 2 else numbers[1]
            # Region number within this record
            region_num = region["qualifiers"]["region_number"][0]
            # antiSMASH-style region name: {record_index}.{region_number}
            region_name = f"{record_index}.{region_num}"
            region_dict = {
                "file": path.stem,
                "record_id": record["name"],
                "record_index": record_index,
                "region": region_num,
                "region_name": region_name,
                "start": start,
                "end": end,
                "contig_edge": region["qualifiers"]["contig_edge"][0],
                "product": " / ".join(region["qualifiers"]["product"]),
                "record_desc": record["description"],
            }
            kcb_dict = {"KCB_hit": "", "KCB_acc": "", "KCB_sim": ""}
            if knownblast:
                hits = knownblast[i]["ranking"]
                if hits:
                    sim = hits[0][1]["similarity"]
                    if sim > KCB_THRESHOLDS['low']:
                        if sim > KCB_THRESHOLDS['high']:
                            sim_level = "high"
                        elif sim > KCB_THRESHOLDS['medium']:
                            sim_level = "medium"
                        else:
                            sim_level = "low"
                        kcb_dict = {
                            "KCB_hit": hits[0][0]["description"],
                            "KCB_acc": hits[0][0]["accession"],
                            "KCB_sim": sim_level,
                        }
            region_dict.update(kcb_dict)
            result_list.append(region_dict)
    return result_list


def main(asdir, outpath):
    record_infos = []

    jsons = asdir.glob("*/*.json")
    for path in jsons:
        record_infos.extend(parse_json(path))

    fieldnames = [
        "file",
        "record_id",
        "record_index",
        "region",
        "region_name",
        "start",
        "end",
        "contig_edge",
        "product",
        "KCB_hit",
        "KCB_acc",
        "KCB_sim",
        "record_desc",
    ]
    with Path.open(outpath, "w") as outf:
        writer = csv.DictWriter(outf, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(record_infos)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Given a bunch of antismash results, tabulate BGC regions"
    )

    parser.add_argument(
        "asdir", type=Path, help="directory containing antiSMASH directories"
    )
    parser.add_argument(
        "outpath", type=Path, help="desired path+name for the output TSV"
    )

    args = parser.parse_args()

    main(args.asdir, args.outpath)
