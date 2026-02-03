#!/usr/bin/env python3
"""Shared antiSMASH JSON parsing utilities."""

import json
import re
from pathlib import Path
from .constants import GENE_COLORS


def load_antismash_json(json_path):
    """Load and return antiSMASH JSON data."""
    try:
        with open(json_path) as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"Warning: Could not read {json_path}: {e}")
        return None


def find_antismash_json(genome_name, antismash_dir):
    """Find antiSMASH JSON file for a genome."""
    antismash_path = Path(antismash_dir)

    # Try exact match
    json_path = antismash_path / genome_name / f"{genome_name}.json"
    if json_path.exists():
        return str(json_path)

    # Try glob for partial match
    for subdir in antismash_path.iterdir():
        if subdir.is_dir() and genome_name in subdir.name:
            json_files = list(subdir.glob("*.json"))
            if json_files:
                return str(json_files[0])

    # Try all subdirectories
    for json_file in antismash_path.glob("*/*.json"):
        if genome_name in json_file.stem or genome_name in str(json_file.parent):
            return str(json_file)

    return None


def parse_location(location_str):
    """Parse antiSMASH location string to start, end, strand."""
    # Handle formats like "[1256:1500](+)", "complement([100:500])", "join{...}"
    numbers = re.findall(r'\d+', location_str)
    if len(numbers) >= 2:
        start = int(min(numbers, key=int))
        end = int(max(numbers, key=int))
    else:
        start = end = int(numbers[0]) if numbers else 0

    strand = -1 if 'complement' in location_str or '(-)' in location_str else 1
    return start, end, strand


def get_gene_function_category(gene_functions, sec_met_domains=None):
    """Determine gene function category from antiSMASH annotations."""
    if not gene_functions:
        return 'other'

    func_str = ' '.join(gene_functions).lower()
    if 'biosynthetic' in func_str:
        if 'core' in func_str or 'rule-based' in func_str:
            return 'biosynthetic'
        return 'biosynthetic-additional'
    if 'regulatory' in func_str:
        return 'regulatory'
    if 'transport' in func_str:
        return 'transport'
    if 'resistance' in func_str:
        return 'resistance'
    return 'other'


def get_gene_color(function_category):
    """Get color for a gene function category."""
    return GENE_COLORS.get(function_category, GENE_COLORS['other'])


def extract_region_info(antismash_data, record_id, region_number):
    """Extract region information from antiSMASH JSON data.

    Returns (region_start, region_end, product) or (None, None, None) if not found.
    """
    if not antismash_data:
        return None, None, None

    # Find the correct record
    target_record = None
    for record in antismash_data.get('records', []):
        if record.get('name') == record_id or record.get('id') == record_id:
            target_record = record
            break

    if not target_record:
        # Try matching by partial name
        for record in antismash_data.get('records', []):
            if record_id in record.get('name', '') or record_id in record.get('id', ''):
                target_record = record
                break

    if not target_record:
        return None, None, None

    # Find the region feature
    for feature in target_record.get('features', []):
        if feature.get('type') == 'region':
            quals = feature.get('qualifiers', {})
            if quals.get('region_number', [''])[0] == str(region_number):
                start, end, _ = parse_location(feature.get('location', ''))
                product = ' / '.join(quals.get('product', ['unknown']))
                return start, end, product

    return None, None, None


def extract_genes_from_region(antismash_data, record_id, region_number):
    """Extract gene/CDS features within a specific BGC region.

    Returns (genes_list, region_start, region_end, product, enhanced_analysis)
    """
    genes = []
    region_start, region_end, region_product = extract_region_info(
        antismash_data, record_id, region_number
    )

    if region_start is None:
        return genes, None, None, None, None

    # Find the correct record
    target_record = None
    for record in antismash_data.get('records', []):
        if record.get('name') == record_id or record.get('id') == record_id:
            target_record = record
            break
        if record_id in record.get('name', '') or record_id in record.get('id', ''):
            target_record = record
            break

    if not target_record:
        return genes, region_start, region_end, region_product, None

    # First pass: collect all locus tags in region
    region_locus_tags = []
    for feature in target_record.get('features', []):
        if feature.get('type') != 'CDS':
            continue
        start, end, strand = parse_location(feature.get('location', ''))
        if start >= region_start and end <= region_end:
            quals = feature.get('qualifiers', {})
            locus_tag = quals.get('locus_tag', ['unknown'])[0]
            region_locus_tags.append(locus_tag)

    # Extract domain hits
    clusterhmmer_by_gene, tigrfam_by_gene = extract_domain_hits(
        target_record, region_locus_tags
    )

    # Extract CDS features within region bounds
    for feature in target_record.get('features', []):
        if feature.get('type') != 'CDS':
            continue

        start, end, strand = parse_location(feature.get('location', ''))

        if start >= region_start and end <= region_end:
            quals = feature.get('qualifiers', {})
            gene_functions = quals.get('gene_functions', [])
            sec_met_domains = quals.get('sec_met_domain', [])
            function_category = get_gene_function_category(gene_functions, sec_met_domains)
            locus_tag = quals.get('locus_tag', ['unknown'])[0]
            translation = quals.get('translation', [''])[0]

            gene_data = {
                'locus_tag': locus_tag,
                'start': start,
                'end': end,
                'strand': strand,
                'product': quals.get('product', ['hypothetical protein'])[0],
                'gene_name': quals.get('gene', [''])[0],
                'function': function_category,
                'gene_functions': gene_functions,
                'domains': sec_met_domains,
                'color': get_gene_color(function_category),
                'clusterhmmer_hits': clusterhmmer_by_gene.get(locus_tag, []),
                'tigrfam_hits': tigrfam_by_gene.get(locus_tag, []),
                'translation': translation
            }
            genes.append(gene_data)

    # Sort genes by start position
    genes.sort(key=lambda x: x['start'])

    # Build enhanced_analysis summary
    all_ch_hits = [hit for hits in clusterhmmer_by_gene.values() for hit in hits]
    all_tf_hits = [hit for hits in tigrfam_by_gene.values() for hit in hits]
    enhanced_analysis = None
    if all_ch_hits or all_tf_hits:
        enhanced_analysis = {
            'clusterhmmer_hits': all_ch_hits,
            'tigrfam_hits': all_tf_hits
        }

    return genes, region_start, region_end, region_product, enhanced_analysis


def extract_domain_hits(record, gene_locus_tags):
    """Extract ClusterHmmer and TIGRFam hits mapped to genes."""
    clusterhmmer_by_gene = {lt: [] for lt in gene_locus_tags}
    tigrfam_by_gene = {lt: [] for lt in gene_locus_tags}

    modules = record.get('modules', {})

    # ClusterHmmer - Pfam domain analysis
    clusterhmmer = modules.get('antismash.detection.cluster_hmmer', {})
    if clusterhmmer:
        for hit in clusterhmmer.get('hits', []):
            locus_tag = hit.get('locus_tag', hit.get('label', ''))
            if locus_tag in clusterhmmer_by_gene:
                clusterhmmer_by_gene[locus_tag].append({
                    'hmm': hit.get('domain', ''),
                    'identifier': hit.get('identifier', ''),
                    'description': hit.get('description', ''),
                    'score': hit.get('score', 0),
                    'evalue': hit.get('evalue', 1),
                    'start': hit.get('protein_start', 0),
                    'end': hit.get('protein_end', 0)
                })

    # TIGRFam analysis
    tigrfam = modules.get('antismash.detection.tigrfam', {})
    if tigrfam:
        tigrfam_hits = tigrfam.get('hits', [])
        if isinstance(tigrfam_hits, list):
            for hit in tigrfam_hits:
                locus_tag = hit.get('locus_tag', hit.get('label', ''))
                if locus_tag in tigrfam_by_gene:
                    tigrfam_by_gene[locus_tag].append({
                        'tigrfam_id': hit.get('identifier', ''),
                        'domain': hit.get('domain', ''),
                        'description': hit.get('description', ''),
                        'score': hit.get('score', 0),
                        'evalue': hit.get('evalue', 1)
                    })

    return clusterhmmer_by_gene, tigrfam_by_gene


def find_record_index(antismash_data, record_id):
    """Find the 1-based index of a record in antiSMASH data for building anchors."""
    if not antismash_data:
        return 1

    for i, record in enumerate(antismash_data.get('records', [])):
        rec_id = record.get('id', '')
        rec_name = record.get('name', '')
        if (record_id in rec_id or record_id in rec_name or
            rec_id.startswith(record_id) or rec_name.startswith(record_id)):
            return i + 1  # 1-based index for antiSMASH anchors
    return 1


def extract_kcb_hits(antismash_data):
    """Extract KnownClusterBlast hits from antiSMASH JSON data.

    Returns list of dicts with keys: record_id, region, hit_name, accession, similarity
    """
    hits = []
    if not antismash_data:
        return hits

    for record in antismash_data.get('records', []):
        record_id = record.get('id', record.get('name', ''))
        modules = record.get('modules', {})
        kcb_module = modules.get('antismash.modules.clusterblast', {})

        if not kcb_module:
            continue

        for region_key, region_data in kcb_module.items():
            if not region_key.startswith('region'):
                continue

            try:
                region_num = int(region_key.replace('region', ''))
            except ValueError:
                continue

            knowncluster = region_data.get('knowncluster', {})
            for ranking in knowncluster.get('rankings', []):
                if not ranking:
                    continue
                hit_info = ranking[0] if isinstance(ranking, list) else ranking
                accession = hit_info.get('accession', '')
                description = hit_info.get('description', '')

                # Calculate similarity from hits
                similarity = 0
                if 'hits' in hit_info:
                    similarity = len(hit_info['hits'])

                hits.append({
                    'record_id': record_id,
                    'region': region_num,
                    'hit_name': description,
                    'accession': accession,
                    'similarity': similarity
                })

    return hits


def build_record_index_map(antismash_dir):
    """Build a mapping from (genome, record_id) -> record_index by scanning antiSMASH JSONs.

    This is useful for constructing region_name identifiers and antiSMASH anchor links.

    Args:
        antismash_dir: Path to directory containing antiSMASH result subdirectories

    Returns:
        dict mapping (genome_name, record_id) -> record_index (1-based)
    """
    record_index_map = {}
    antismash_path = Path(antismash_dir)

    for json_file in antismash_path.glob("*/*.json"):
        genome_name = json_file.parent.name
        try:
            with open(json_file) as f:
                data = json.load(f)
            for i, record in enumerate(data.get('records', []), 1):
                rec_id = record.get('id', '')
                rec_name = record.get('name', '')
                # Store both id and name as keys (some may have .1 suffix)
                if rec_id:
                    record_index_map[(genome_name, rec_id)] = i
                    # Also store without .1 suffix for matching
                    if rec_id.endswith('.1'):
                        record_index_map[(genome_name, rec_id[:-2])] = i
                if rec_name and rec_name != rec_id:
                    record_index_map[(genome_name, rec_name)] = i
        except Exception:
            continue

    return record_index_map
