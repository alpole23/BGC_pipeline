#!/usr/bin/env python3
"""
Build phylogenetic trees for phosphonate BGC coupling enzyme analysis.

Two tree types:
  A: Combined pepM tree — universal marker (SMCOG1231 / PEP_mutase), all BGCs + reference
     anchors. Annotated with coupling enzyme class, GCF family, and source (ref vs query).
     Unknown BGCs fall naturally into the closest clade.

  B: Per-class coupling enzyme trees — one tree per class using the class-defining marker
     gene. Ppd and Ppd-CDP share one tree (same enzyme; class distinction is an annotation
     layer). References anchor each class tree.

HMM strategy (no external MSA tool required):
  1. hmmbuild from single seed reference    → initial HMM
  2. hmmalign all references to initial HMM → aligned references
  3. hmmbuild from aligned references       → refined HMM
  4. hmmalign all (refs + queries)          → final alignment
  5. BioPython NJ (BLOSUM62)               → .nwk tree

Usage:
    python scripts/bgc_coupling_tree.py \\
        --antismash_dir  results/antismash_results/Pantoea \\
        --metadata       results/bgc_trees/Pantoea/phosphonate_metadata.json \\
        --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \\
        --ref_pepm_faa   results/bgc_trees/Pantoea/coupling_enzyme_trees/reference_pepM.faa \\
        --ref_coupling_faa results/bgc_trees/Pantoea/coupling_enzyme_trees/reference_coupling_enzymes.faa \\
        --outdir         results/bgc_trees/Pantoea/coupling_enzyme_trees \\
        --hmmbuild       /path/to/hmmbuild \\
        --hmmalign       /path/to/hmmalign \\
        --tree           both
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent))
from utils.constants import COUPLING_COLORS as _BASE_COUPLING_COLORS, LEGACY_CLASS_NAMES

# Add the Reference class used only in coupling trees
COUPLING_COLORS = {**_BASE_COUPLING_COLORS, 'Reference': '#333333'}

# Coupling enzyme class for each reference pepM BGC
REF_PEPM_CLASS = {
    'BGC0000904':          'FrbC-like',   # FR-900098, FrbD pepM
    'BGC0000897':          'Ppd',         # Dehydrophos, DhpE pepM
    'BGC0000938':          'Ppd',         # Fosfomycin, Fom1 pepM
    'BGC0000806':          'Ppd',         # 2-AEP, Glycomyces pepM
    'Phosphonoalamide_BGC':'PalB-like',   # PnaD pepM
    'Valinophos_BGC':      'VlpB-like',   # VlpA pepM
    'Pantaphos_BGC':       'FrbC-like',   # Pantaphos, HvrA pepM
}

# Coupling enzyme class for each reference coupling enzyme entry (by protein name)
REF_COUPLING_CLASS = {
    'FrbC': 'FrbC-like',
    'HvrC': 'FrbC-like',
    'DhpF': 'Ppd',
    'Fom2': 'Ppd',
    'Ppd':  'Ppd',
    'VlpB': 'VlpB-like',
    'PnaA': 'PalB-like',
}

# SMCOG / domain markers for coupling enzyme CDS extraction
CLASS_MARKERS = {
    'FrbC-like': ('smcog',  'SMCOG1271'),
    'Ppd':       ('smcog',  'SMCOG1055'),
    'Ppd-CDP':   ('smcog',  'SMCOG1055'),
    'VlpB-like': ('domain', 'Fe-ADH'),
    'PalB-like': ('smcog',  'SMCOG1013'),
}

GCF_PALETTE = [
    '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e',
    '#e6ab02', '#a6761d', '#666666', '#a6cee3', '#b2df8a',
]


# ─── Data loading ─────────────────────────────────────────────────────────────


def load_coupling_classes(path):
    """Read iTOL colorstrip file → {bgc_label: class_id}."""
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
                    cls = parts[2]
                    cls = LEGACY_CLASS_NAMES.get(cls, cls)
                    classes[parts[0]] = cls
    return classes


def load_metadata(path):
    """Read BGC metadata JSON → {label: {organism, gcf, gbk_path}}."""
    with open(path) as f:
        records = json.load(f)
    meta = {}
    for r in records:
        lbl = r['label']
        gcf = next(
            (f['family_id'] for f in r.get('families', []) if f['cutoff'] == 0.3),
            None,
        )
        meta[lbl] = {
            'organism': r.get('organism', ''),
            'gcf':      gcf,
            'gbk_path': r.get('gbk_path', ''),
        }
    return meta


def build_json_index(antismash_dir):
    """Walk antismash_dir → {genome_name: json_path}."""
    index = {}
    for genome in os.listdir(antismash_dir):
        jp = os.path.join(antismash_dir, genome, f'{genome}.json')
        if os.path.exists(jp):
            index[genome] = jp
    return index


def genome_from_gbk_path(gbk_path):
    return os.path.basename(os.path.dirname(gbk_path))


def parse_label(label):
    """'CONTIG.regionNNN' → (contig_id, zero-padded region str)."""
    label = re.sub(r'_\d+$', '', label)
    m = re.search(r'^(.+?)\.region(\d+)$', label)
    if m:
        return m.group(1), m.group(2).zfill(3)
    return label, '001'


def parse_location_bounds(loc_str):
    """Extract (min_start, max_end) from an antiSMASH location string."""
    coords = re.findall(r'\[(\d+):(\d+)\]', str(loc_str))
    if not coords:
        return None, None
    starts = [int(s) for s, _ in coords]
    ends   = [int(e) for _, e in coords]
    return min(starts), max(ends)


# ─── Sequence extraction ──────────────────────────────────────────────────────

def _cds_in_region(feat, region_start, region_end):
    """Return True if a CDS feature overlaps the region."""
    if region_start is None:
        return True
    start, end = parse_location_bounds(feat.get('location', ''))
    if start is None:
        return True
    return not (end < region_start or start > region_end)


def _get_region_bounds(rec, contig_id, region_num):
    """Find the phosphonate region and return its (start, end)."""
    for feat in rec.get('features', []):
        if feat.get('type') != 'region':
            continue
        rnum = str(feat.get('qualifiers', {}).get('region_number', ['?'])[0]).zfill(3)
        products = feat.get('qualifiers', {}).get('product', [])
        if rnum == region_num and any('phosphonate' in p for p in products):
            return parse_location_bounds(feat.get('location', ''))
    return None, None


def extract_cds_from_json(json_path, contig_id, region_num,
                           is_pepm=False, smcog=None, domain=None):
    """
    Find and return (translation, gene_name) for a CDS within the specified
    phosphonate region matching the given marker.

    Specify one of:
      is_pepm=True       — match SMCOG1231 or PEP_mutase sec_met_domain
      smcog='SMCOG1271'  — match a specific SMCOG annotation
      domain='Fe-ADH'    — match a specific sec_met_domain or rule-based-cluster
    """
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception:
        return None, None

    for rec in data['records']:
        if contig_id not in rec.get('id', ''):
            continue

        r_start, r_end = _get_region_bounds(rec, contig_id, region_num)

        for feat in rec.get('features', []):
            if feat.get('type') != 'CDS':
                continue
            if not _cds_in_region(feat, r_start, r_end):
                continue

            quals = feat.get('qualifiers', {})
            translation = quals.get('translation', [''])[0]
            if not translation:
                continue

            gene_functions = quals.get('gene_functions', [])
            sec_met        = quals.get('sec_met_domain', [])

            matched = False
            if is_pepm:
                matched = (
                    any('SMCOG1231' in gf for gf in gene_functions) or
                    any('PEP_mutase' in sd for sd in sec_met)
                )
            elif smcog:
                matched = any(smcog in gf for gf in gene_functions)
            elif domain:
                matched = (
                    any(domain in sd for sd in sec_met) or
                    any('rule-based-clusters' in gf and domain in gf
                        for gf in gene_functions)
                )

            if matched:
                gene = quals.get('gene', [''])[0] or quals.get('locus_tag', [''])[0]
                return translation, gene

    return None, None


# ─── HMM build and alignment ──────────────────────────────────────────────────

def write_fasta(records, path):
    with open(path, 'w') as f:
        SeqIO.write(records, f, 'fasta')


def run_hmmbuild(input_faa, out_hmm, hmmbuild_bin, name='profile'):
    cmd = [hmmbuild_bin, '--amino', '-n', name, out_hmm, input_faa]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  hmmbuild failed:\n{r.stderr[:600]}', file=sys.stderr)
        return False
    return True


def run_hmmalign(hmm, seqs_faa, out_afa, hmmalign_bin):
    cmd = [hmmalign_bin, '--amino', '--trim', '--outformat', 'afa', hmm, seqs_faa]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f'  hmmalign failed:\n{r.stderr[:600]}', file=sys.stderr)
        return False
    with open(out_afa, 'w') as f:
        f.write(r.stdout)
    return True


def hmm_align_all(ref_records, query_records, workdir, hmmbuild_bin, hmmalign_bin, label):
    """
    Build a refined HMM from references and align refs + queries to it.
    Returns (MultipleSeqAlignment, aligned_faa_path) or (None, None) on failure.
    """
    os.makedirs(workdir, exist_ok=True)

    seed_faa  = os.path.join(workdir, 'seed.faa')
    refs_faa  = os.path.join(workdir, 'refs.faa')
    all_faa   = os.path.join(workdir, 'all.faa')
    hmm_init  = os.path.join(workdir, 'init.hmm')
    refs_afa  = os.path.join(workdir, 'refs_aligned.faa')
    hmm_final = os.path.join(workdir, 'final.hmm')
    all_afa   = os.path.join(workdir, 'all_aligned.faa')

    write_fasta([ref_records[0]], seed_faa)
    write_fasta(ref_records, refs_faa)
    write_fasta(ref_records + query_records, all_faa)

    print(f'    [1/4] hmmbuild from seed: {ref_records[0].id}')
    if not run_hmmbuild(seed_faa, hmm_init, hmmbuild_bin, label + '_init'):
        return None, None

    print(f'    [2/4] hmmalign {len(ref_records)} references...')
    if not run_hmmalign(hmm_init, refs_faa, refs_afa, hmmalign_bin):
        return None, None

    print(f'    [3/4] hmmbuild refined HMM from aligned references...')
    if not run_hmmbuild(refs_afa, hmm_final, hmmbuild_bin, label):
        return None, None

    print(f'    [4/4] hmmalign {len(ref_records) + len(query_records)} sequences...')
    if not run_hmmalign(hmm_final, all_faa, all_afa, hmmalign_bin):
        return None, None

    alignment = AlignIO.read(all_afa, 'fasta')
    print(f'    Alignment: {len(alignment)} seqs × {alignment.get_alignment_length()} cols')
    return alignment, all_afa


# ─── Tree building ─────────────────────────────────────────────────────────────

def build_nj_tree(alignment, out_nwk):
    # hmmalign uses '.' for insert-state gaps, '-' for match-state deletions,
    # and lowercase letters for residues aligned to insert states.
    # BioPython DistanceCalculator requires uppercase residues and '-' only.
    from Bio.Align import MultipleSeqAlignment
    cleaned = MultipleSeqAlignment([
        SeqRecord(Seq(str(r.seq).upper().replace('.', '-')), id=r.id, description='')
        for r in alignment
    ])
    print(f'    Computing BLOSUM62 distance matrix ({len(cleaned)} × {len(cleaned)})...')
    dm = DistanceCalculator('blosum62').get_distance(cleaned)
    print(f'    Building NJ tree...')
    tree = DistanceTreeConstructor().nj(dm)
    # BioPython's Newick writer is recursive; large trees exceed the default limit.
    import sys
    old_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old_limit, len(alignment) * 10))
    try:
        with open(out_nwk, 'w') as f:
            Phylo.write(tree, f, 'newick')
    finally:
        sys.setrecursionlimit(old_limit)
    print(f'    Written: {out_nwk}')
    return tree


# ─── iTOL annotation writers ──────────────────────────────────────────────────

def write_itol_colorstrip(labels, color_fn, display_fn, dataset_label, legend_items, out_path):
    """
    Generic DATASET_COLORSTRIP writer.
    color_fn(label)   → hex color string
    display_fn(label) → display label string (shown in strip tooltip)
    legend_items      → [(label_str, color_str), ...]
    """
    with open(out_path, 'w') as f:
        f.write('DATASET_COLORSTRIP\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{dataset_label}\n')
        f.write('COLOR\t#333333\n')
        f.write('STRIP_WIDTH\t40\n')
        f.write('SHOW_BORDER\t1\n')
        f.write('BORDER_WIDTH\t0.5\n')
        if legend_items:
            f.write(f'LEGEND_TITLE\t{dataset_label}\n')
            f.write('LEGEND_SHAPES\t' + '\t'.join('1' for _ in legend_items) + '\n')
            f.write('LEGEND_COLORS\t' + '\t'.join(c for _, c in legend_items) + '\n')
            f.write('LEGEND_LABELS\t' + '\t'.join(l for l, _ in legend_items) + '\n')
        f.write('DATA\n')
        for lbl in labels:
            f.write(f'{lbl}\t{color_fn(lbl)}\t{display_fn(lbl)}\n')


def write_itol_text(labels, text_fn, dataset_label, out_path):
    """DATASET_TEXT writer for leaf label annotations (e.g. organism names)."""
    with open(out_path, 'w') as f:
        f.write('DATASET_TEXT\n')
        f.write('SEPARATOR TAB\n')
        f.write(f'DATASET_LABEL\t{dataset_label}\n')
        f.write('COLOR\t#333333\n')
        f.write('DATA\n')
        for lbl in labels:
            text = text_fn(lbl)
            if text:
                # node_id, text, position (1=after), color, style, size_factor
                f.write(f'{lbl}\t{text}\t1\t#333333\tnormal\t1\n')


def write_tree_a_itol(seq_labels, coupling_classes, metadata, ref_records, outdir):
    """Write all iTOL annotation files for Tree A."""

    # Determine coupling class for each label (including references)
    def get_class(lbl):
        if lbl.startswith('REF|'):
            bgc_id = lbl.split('|')[1]  # e.g. BGC0000904
            return REF_PEPM_CLASS.get(bgc_id, 'Reference')
        return coupling_classes.get(lbl, 'Unknown')

    present_classes = sorted(set(get_class(l) for l in seq_labels))
    coupling_legend = [(cls, COUPLING_COLORS.get(cls, '#aaaaaa')) for cls in present_classes]

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: COUPLING_COLORS.get(get_class(l), '#aaaaaa'),
        display_fn = get_class,
        dataset_label = 'Coupling enzyme class',
        legend_items  = coupling_legend,
        out_path = os.path.join(outdir, 'itol_coupling_class.txt'),
    )

    # GCF colorstrip
    gcf_set = sorted(set(
        metadata[l]['gcf'] for l in seq_labels
        if l in metadata and metadata[l]['gcf'] is not None
    ))
    gcf_color = {g: GCF_PALETTE[i % len(GCF_PALETTE)] for i, g in enumerate(gcf_set)}

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: gcf_color.get(metadata[l]['gcf'], '#dddddd')
                               if l in metadata and metadata[l]['gcf'] is not None
                               else '#333333' if l.startswith('REF|') else '#dddddd',
        display_fn = lambda l: f'GCF-{metadata[l]["gcf"]}'
                               if l in metadata and metadata[l]['gcf'] is not None
                               else 'Reference' if l.startswith('REF|') else 'No GCF',
        dataset_label = 'GCF family',
        legend_items  = [(f'GCF-{g}', gcf_color[g]) for g in gcf_set],
        out_path = os.path.join(outdir, 'itol_gcf.txt'),
    )

    # Source (reference vs query) colorstrip
    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: '#333333' if l.startswith('REF|') else '#cccccc',
        display_fn = lambda l: 'Reference' if l.startswith('REF|') else 'Query',
        dataset_label = 'Source',
        legend_items  = [('Reference', '#333333'), ('Query (Pantoea)', '#cccccc')],
        out_path = os.path.join(outdir, 'itol_source.txt'),
    )

    # Text labels — organism name
    ref_map = {}
    for r in ref_records:
        parts = r.id.split('|')  # BGC|acc|name|function|organism
        org = parts[4].replace('_', ' ') if len(parts) > 4 else r.id
        ref_map[f'REF|{r.id}'] = org

    write_itol_text(
        seq_labels,
        text_fn = lambda l: ref_map.get(l) or
                            metadata.get(l, {}).get('organism', ''),
        dataset_label = 'Organism',
        out_path = os.path.join(outdir, 'itol_organism.txt'),
    )


def write_tree_b_itol(seq_labels, coupling_classes, metadata, ref_records, outdir):
    """Write iTOL annotation files for a Tree B class subtree."""

    # For Tree B, coupling class is the class being analyzed; refs are anchors
    def get_class(lbl):
        if lbl.startswith('REF|'):
            return 'Reference'
        return coupling_classes.get(lbl, 'Unknown')

    present = sorted(set(get_class(l) for l in seq_labels))
    coupling_legend = [(cls, COUPLING_COLORS.get(cls, '#aaaaaa')) for cls in present]

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: COUPLING_COLORS.get(get_class(l), '#aaaaaa'),
        display_fn = get_class,
        dataset_label = 'Coupling enzyme class',
        legend_items  = coupling_legend,
        out_path = os.path.join(outdir, 'itol_coupling_class.txt'),
    )

    gcf_set = sorted(set(
        metadata[l]['gcf'] for l in seq_labels
        if l in metadata and metadata[l]['gcf'] is not None
    ))
    gcf_color = {g: GCF_PALETTE[i % len(GCF_PALETTE)] for i, g in enumerate(gcf_set)}

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: gcf_color.get(metadata[l]['gcf'], '#dddddd')
                               if l in metadata and metadata[l]['gcf'] is not None
                               else '#333333' if l.startswith('REF|') else '#dddddd',
        display_fn = lambda l: f'GCF-{metadata[l]["gcf"]}'
                               if l in metadata and metadata[l]['gcf'] is not None
                               else 'Reference' if l.startswith('REF|') else 'No GCF',
        dataset_label = 'GCF family',
        legend_items  = [(f'GCF-{g}', gcf_color[g]) for g in gcf_set],
        out_path = os.path.join(outdir, 'itol_gcf.txt'),
    )

    write_itol_colorstrip(
        seq_labels,
        color_fn   = lambda l: '#333333' if l.startswith('REF|') else '#cccccc',
        display_fn = lambda l: 'Reference' if l.startswith('REF|') else 'Query',
        dataset_label = 'Source',
        legend_items  = [('Reference', '#333333'), ('Query (Pantoea)', '#cccccc')],
        out_path = os.path.join(outdir, 'itol_source.txt'),
    )

    ref_map = {}
    for r in ref_records:
        parts = r.id.split('|')
        org = parts[4].replace('_', ' ') if len(parts) > 4 else r.id
        ref_map[f'REF|{r.id}'] = org

    write_itol_text(
        seq_labels,
        text_fn = lambda l: ref_map.get(l) or metadata.get(l, {}).get('organism', ''),
        dataset_label = 'Organism',
        out_path = os.path.join(outdir, 'itol_organism.txt'),
    )


# ─── Tree A ────────────────────────────────────────────────────────────────────

def build_tree_a(args, metadata, coupling_classes, json_index, ref_pepm_records):
    print('\n=== Tree A: Combined pepM tree ===')
    outdir = os.path.join(args.outdir, 'tree_A')
    os.makedirs(outdir, exist_ok=True)

    print(f'Extracting pepM sequences from {len(metadata)} BGCs...')
    query_records = []
    missing = []

    for lbl, meta in metadata.items():
        contig_id, region_num = parse_label(lbl)
        genome    = genome_from_gbk_path(meta['gbk_path']) if meta['gbk_path'] else None
        json_path = json_index.get(genome)

        seq, _ = extract_cds_from_json(json_path, contig_id, region_num, is_pepm=True) \
                 if json_path else (None, None)

        if seq:
            query_records.append(SeqRecord(Seq(seq), id=lbl, description=''))
        else:
            missing.append(lbl)

    print(f'  Extracted: {len(query_records)}  |  Missing pepM: {len(missing)}')
    if missing:
        missing_log = os.path.join(outdir, 'missing_pepm.txt')
        with open(missing_log, 'w') as f:
            f.write('\n'.join(missing))
        print(f'  Missing labels written to: {missing_log}')

    ref_records = [
        SeqRecord(Seq(str(r.seq)), id=f'REF|{r.id}', description=r.description)
        for r in ref_pepm_records
    ]

    alignment, _ = hmm_align_all(
        ref_records, query_records,
        os.path.join(outdir, 'hmm_work'),
        args.hmmbuild, args.hmmalign, 'pepm',
    )
    if alignment is None:
        print('Tree A: alignment failed.', file=sys.stderr)
        return

    seq_labels = [r.id for r in alignment]
    build_nj_tree(alignment, os.path.join(outdir, 'pepm_tree.nwk'))

    print('  Writing iTOL annotations...')
    write_tree_a_itol(seq_labels, coupling_classes, metadata, ref_pepm_records, outdir)
    print(f'Tree A complete → {outdir}')


# ─── Tree B ────────────────────────────────────────────────────────────────────

def build_tree_b(args, metadata, coupling_classes, json_index, ref_coupling_records):
    print('\n=== Tree B: Per-class coupling enzyme trees ===')
    outdir_b = os.path.join(args.outdir, 'tree_B')
    os.makedirs(outdir_b, exist_ok=True)

    # Group BGCs by class; merge Ppd + Ppd-CDP into one tree
    class_bgcs = defaultdict(list)
    for lbl, cls in coupling_classes.items():
        key = 'Ppd+Ppd-CDP' if cls in ('Ppd', 'Ppd-CDP') else cls
        class_bgcs[key].append(lbl)

    # Group reference sequences by class; merge Ppd + Ppd-CDP refs
    ref_by_class = defaultdict(list)
    for r in ref_coupling_records:
        parts = r.id.split('|')
        name = parts[2] if len(parts) > 2 else ''
        cls  = REF_COUPLING_CLASS.get(name)
        if cls:
            key = 'Ppd+Ppd-CDP' if cls in ('Ppd', 'Ppd-CDP') else cls
            ref_by_class[key].append(r)

    for class_key, bgc_labels in sorted(class_bgcs.items()):
        if class_key == 'Unknown':
            print(f'\n  Skipping Unknown class (no coupling enzyme to extract)')
            continue

        print(f'\n  Class: {class_key} ({len(bgc_labels)} BGCs)')
        class_outdir = os.path.join(outdir_b, class_key.replace('+', '_'))
        os.makedirs(class_outdir, exist_ok=True)

        # Determine extraction marker
        if class_key == 'Ppd+Ppd-CDP':
            marker_type, marker_value = 'smcog', 'SMCOG1055'
        else:
            marker_type, marker_value = CLASS_MARKERS[class_key]

        # Extract coupling enzyme sequences
        query_records = []
        missing = []
        for lbl in bgc_labels:
            meta      = metadata.get(lbl, {})
            contig_id, region_num = parse_label(lbl)
            genome    = genome_from_gbk_path(meta.get('gbk_path', '')) if meta.get('gbk_path') else None
            json_path = json_index.get(genome)

            seq, _ = extract_cds_from_json(
                json_path, contig_id, region_num,
                smcog=marker_value  if marker_type == 'smcog'  else None,
                domain=marker_value if marker_type == 'domain' else None,
            ) if json_path else (None, None)

            if seq:
                query_records.append(SeqRecord(Seq(seq), id=lbl, description=''))
            else:
                missing.append(lbl)

        print(f'    Extracted: {len(query_records)}  |  Missing: {len(missing)}')

        refs = ref_by_class.get(class_key, [])
        if not refs:
            print(f'    No references for {class_key} — skipping.')
            continue

        ref_records = [
            SeqRecord(Seq(str(r.seq)), id=f'REF|{r.id}', description=r.description)
            for r in refs
        ]

        alignment, _ = hmm_align_all(
            ref_records, query_records,
            os.path.join(class_outdir, 'hmm_work'),
            args.hmmbuild, args.hmmalign, class_key.lower().replace('+', '_'),
        )
        if alignment is None:
            print(f'    Tree B/{class_key}: alignment failed.', file=sys.stderr)
            continue

        seq_labels = [r.id for r in alignment]
        out_nwk = os.path.join(class_outdir, f'{class_key.replace("+", "_")}_tree.nwk')
        build_nj_tree(alignment, out_nwk)

        print(f'    Writing iTOL annotations...')
        write_tree_b_itol(seq_labels, coupling_classes, metadata, refs, class_outdir)
        print(f'    Tree B/{class_key} complete → {class_outdir}')


# ─── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Build pepM and coupling enzyme phylogenetic trees for phosphonate BGCs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--antismash_dir',       required=True,
                        help='antiSMASH results directory')
    parser.add_argument('--metadata',            required=True,
                        help='BGC metadata JSON from bgc_pfam_tree.py')
    parser.add_argument('--coupling_annotation', required=True,
                        help='iTOL colorstrip from bgc_coupling_annotation.py')
    parser.add_argument('--ref_pepm_faa',        required=True,
                        help='Reference pepM sequences (FASTA)')
    parser.add_argument('--ref_coupling_faa',    required=True,
                        help='Reference coupling enzyme sequences (FASTA)')
    parser.add_argument('--outdir',              required=True,
                        help='Output directory')
    parser.add_argument('--hmmbuild',            required=True,
                        help='Path to hmmbuild binary')
    parser.add_argument('--hmmalign',            required=True,
                        help='Path to hmmalign binary')
    parser.add_argument('--tree', default='both', choices=['A', 'B', 'both'],
                        help='Which tree(s) to build (default: both)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print('Loading metadata and annotations...')
    metadata        = load_metadata(args.metadata)
    coupling_classes = load_coupling_classes(args.coupling_annotation)
    print(f'  {len(metadata)} BGCs in metadata')
    print(f'  {len(coupling_classes)} coupling class assignments')

    print('Building antiSMASH JSON index...')
    json_index = build_json_index(args.antismash_dir)
    print(f'  {len(json_index)} genome JSON files')

    ref_pepm_records     = list(SeqIO.parse(args.ref_pepm_faa,     'fasta'))
    ref_coupling_records = list(SeqIO.parse(args.ref_coupling_faa, 'fasta'))
    print(f'  {len(ref_pepm_records)} pepM reference sequences')
    print(f'  {len(ref_coupling_records)} coupling enzyme reference sequences')

    if args.tree in ('A', 'both'):
        build_tree_a(args, metadata, coupling_classes, json_index, ref_pepm_records)

    if args.tree in ('B', 'both'):
        build_tree_b(args, metadata, coupling_classes, json_index, ref_coupling_records)

    print('\nDone.')


if __name__ == '__main__':
    main()
