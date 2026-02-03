"""Shared utilities for BGC analysis pipeline."""

from .constants import BGC_COLORS, GENE_COLORS, GCF_COLORS, KCB_THRESHOLDS, KCB_COLORS
from .parsers import (
    parse_trace_file, parse_duration, parse_memory, format_bytes,
    format_duration_str, parse_timestamp, parse_newick, sanitize_taxon
)
from .antismash_parser import (
    load_antismash_json, find_antismash_json, parse_location,
    get_gene_function_category, get_gene_color, extract_region_info,
    extract_genes_from_region, extract_domain_hits, find_record_index,
    extract_kcb_hits
)
