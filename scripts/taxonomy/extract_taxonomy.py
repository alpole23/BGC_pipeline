#!/usr/bin/env python3
"""Extract taxonomy information using hybrid taxonkit + Entrez approach."""

import argparse
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def get_lineage_taxonkit(tax_ids, taxdump_dir):
    """
    Query taxonkit for lineage information.
    Returns dict mapping taxId -> lineage info
    """
    if not tax_ids:
        return {}

    results = {}
    tax_id_list = [str(tid) for tid in tax_ids]

    try:
        # Run taxonkit lineage with full lineage info
        # Format: taxid, name, lineage, lineage_taxids, lineage_ranks
        proc = subprocess.run(
            ['taxonkit', 'lineage', '-t', '-R',
             '--data-dir', taxdump_dir],
            input='\n'.join(tax_id_list),
            capture_output=True,
            text=True,
            timeout=300
        )

        if proc.returncode != 0:
            print(f"taxonkit warning: {proc.stderr}")

        for line in proc.stdout.strip().split('\n'):
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 4:
                continue

            # taxonkit lineage -t -R output format:
            # taxid \t lineage_names \t lineage_taxids \t lineage_ranks
            tax_id = int(parts[0])
            lineage_names_str = parts[1]  # semicolon-separated names
            lineage_taxids_str = parts[2]  # semicolon-separated taxids
            lineage_ranks_str = parts[3]   # semicolon-separated ranks

            if not lineage_names_str:
                continue

            # Parse lineage into structured format
            lineage_names = lineage_names_str.split(';')
            lineage_ids = lineage_taxids_str.split(';') if lineage_taxids_str else []
            ranks = lineage_ranks_str.split(';') if lineage_ranks_str else []

            # Build structured lineage dict
            lineage = {}
            standard_ranks = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

            for i, rank in enumerate(ranks):
                rank_lower = rank.lower().strip()
                # Map superkingdom to domain (for compatibility)
                if rank_lower == 'superkingdom':
                    rank_lower = 'domain'

                if rank_lower in standard_ranks:
                    lineage[rank_lower] = {
                        'name': lineage_names[i].strip() if i < len(lineage_names) else '',
                        'id': int(lineage_ids[i]) if i < len(lineage_ids) and lineage_ids[i].isdigit() else 0
                    }

            # Scientific name is the last item in the lineage
            scientific_name = lineage_names[-1].strip() if lineage_names else ''

            results[tax_id] = {
                'taxId': tax_id,
                'scientificName': scientific_name,
                'lineage': lineage,
                'rank': ranks[-1].lower().strip() if ranks else 'unknown',
                'source': 'taxonkit'
            }

    except subprocess.TimeoutExpired:
        print("taxonkit timed out")
    except Exception as e:
        print(f"taxonkit error: {e}")

    return results


def get_lineage_entrez(tax_ids, email="pipeline@example.com"):
    """
    Fallback: Query NCBI Entrez for lineage information.
    Handles missing taxIds from taxonkit (new taxa, merged IDs, etc.)
    """
    if not tax_ids:
        return {}

    results = {}

    try:
        from Bio import Entrez
        Entrez.email = email

        # Batch query (max 10000 per request)
        batch_size = 500  # Conservative batch size
        tax_id_list = list(tax_ids)

        for i in range(0, len(tax_id_list), batch_size):
            batch = tax_id_list[i:i+batch_size]
            print(f"  Entrez batch {i//batch_size + 1}: fetching {len(batch)} taxIds...")

            try:
                handle = Entrez.efetch(db="taxonomy", id=batch, retmode="xml")
                records = Entrez.read(handle)
                handle.close()

                for record in records:
                    tax_id = int(record['TaxId'])
                    name = record.get('ScientificName', '')
                    rank = record.get('Rank', 'unknown').lower()

                    # Parse LineageEx for structured lineage
                    lineage = {}
                    for item in record.get('LineageEx', []):
                        item_rank = item.get('Rank', '').lower()
                        if item_rank == 'superkingdom':
                            item_rank = 'domain'

                        if item_rank in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                            lineage[item_rank] = {
                                'name': item.get('ScientificName', ''),
                                'id': int(item.get('TaxId', 0))
                            }

                    # Add the taxon itself if it's at a standard rank
                    if rank in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                        lineage[rank] = {
                            'name': name,
                            'id': tax_id
                        }

                    results[tax_id] = {
                        'taxId': tax_id,
                        'scientificName': name,
                        'lineage': lineage,
                        'rank': rank,
                        'source': 'entrez'
                    }

                # Rate limiting: NCBI allows 3 requests/sec without API key
                time.sleep(0.35)

            except Exception as e:
                print(f"  Entrez batch error: {e}")
                continue

    except ImportError:
        print("Warning: Biopython not available for Entrez fallback")
    except Exception as e:
        print(f"Entrez error: {e}")

    return results


def main(assembly_report, taxdump_dir, output_file):
    print("=" * 60)
    print("TAXONOMY EXTRACTION - Hybrid taxonkit + Entrez approach")
    print("=" * 60)

    # Read assembly report to get all taxIds
    print("\nReading assembly report...")
    assemblies = {}
    tax_ids = set()

    with open(assembly_report, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)

            accession = record.get('accession')
            organism_info = record.get('organism', {})
            tax_id = organism_info.get('taxId')
            organism_name = organism_info.get('organismName', '')

            if not accession or not tax_id:
                continue

            # Extract strain information
            infraspecific = organism_info.get('infraspecificNames', {})
            strain = infraspecific.get('strain', '')

            assemblies[accession] = {
                'taxId': tax_id,
                'organism': organism_name,
                'strain': strain
            }
            tax_ids.add(tax_id)

    print(f"Found {len(assemblies)} assemblies with {len(tax_ids)} unique taxIds")

    # Step 1: Query taxonkit for all taxIds (fast, local)
    print(f"\nStep 1: Querying taxonkit (local database)...")
    taxonkit_results = get_lineage_taxonkit(tax_ids, taxdump_dir)
    print(f"  taxonkit resolved {len(taxonkit_results)} of {len(tax_ids)} taxIds")

    # Step 2: Find missing taxIds and query Entrez (fallback)
    missing_tax_ids = tax_ids - set(taxonkit_results.keys())

    entrez_results = {}
    if missing_tax_ids:
        print(f"\nStep 2: Querying Entrez for {len(missing_tax_ids)} missing taxIds...")
        entrez_results = get_lineage_entrez(missing_tax_ids)
        print(f"  Entrez resolved {len(entrez_results)} additional taxIds")

    # Combine results
    taxonomy_db = {**taxonkit_results, **entrez_results}

    # Report any still-unresolved taxIds
    still_missing = tax_ids - set(taxonomy_db.keys())
    if still_missing:
        print(f"\nWarning: {len(still_missing)} taxIds could not be resolved:")
        for tid in list(still_missing)[:10]:
            print(f"  - {tid}")
        if len(still_missing) > 10:
            print(f"  ... and {len(still_missing) - 10} more")

    # Build taxonomy map for all assemblies
    print("\nBuilding taxonomy map...")
    taxonomy_map = {}

    for accession, assembly_info in assemblies.items():
        tax_id = assembly_info['taxId']
        tax_entry = taxonomy_db.get(tax_id, {})

        taxonomy_map[accession] = {
            'assembly_id': accession,
            'taxId': tax_id,
            'organism': assembly_info['organism'],
            'strain': assembly_info['strain'],
            'lineage': tax_entry.get('lineage', {}),
            'rank': tax_entry.get('rank', 'unknown'),
            'scientificName': tax_entry.get('scientificName', assembly_info['organism']),
            'taxonomy_source': tax_entry.get('source', 'unresolved')
        }

    # Write taxonomy map
    print("\nWriting taxonomy_map.json...")
    with open(output_file, 'w') as f:
        json.dump(taxonomy_map, f, indent=2)

    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total assemblies: {len(taxonomy_map)}")

    source_counts = defaultdict(int)
    rank_counts = defaultdict(int)
    lineage_completeness = defaultdict(int)

    for entry in taxonomy_map.values():
        source_counts[entry.get('taxonomy_source', 'unknown')] += 1
        rank_counts[entry.get('rank', 'unknown')] += 1
        lineage = entry.get('lineage', {})
        for level in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            if level in lineage:
                lineage_completeness[level] += 1

    print("\nTaxonomy source:")
    for source, count in sorted(source_counts.items()):
        print(f"  {source}: {count}")

    print("\nLineage completeness:")
    for level in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        count = lineage_completeness.get(level, 0)
        pct = count / len(taxonomy_map) * 100 if taxonomy_map else 0
        print(f"  {level}: {count} ({pct:.1f}%)")

    print("\nDone!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract taxonomy using hybrid taxonkit + Entrez approach")
    parser.add_argument("assembly_report", help="Path to assembly data report (JSONL)")
    parser.add_argument("taxdump_dir", help="Path to taxonkit taxdump directory")
    parser.add_argument("output_file", help="Path to output JSON file")
    args = parser.parse_args()

    main(args.assembly_report, args.taxdump_dir, args.output_file)
