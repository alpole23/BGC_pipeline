#!/usr/bin/env python3
"""Shared parsing utilities for BGC analysis pipeline."""

import re
from pathlib import Path
from datetime import datetime


def parse_trace_file(trace_path):
    """Parse Nextflow trace TSV file."""
    tasks = []
    if not trace_path or not Path(trace_path).exists():
        return tasks

    with open(trace_path, 'r') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                header = parts
                continue
            task = dict(zip(header, parts))
            tasks.append(task)
    return tasks


def parse_duration(duration_str):
    """Parse Nextflow duration string to seconds."""
    if not duration_str or duration_str == '-':
        return 0
    total_seconds = 0
    if 'ms' in duration_str:
        match = re.search(r'([\d.]+)ms', duration_str)
        if match:
            return float(match.group(1)) / 1000
    match = re.search(r'(\d+)h', duration_str)
    if match:
        total_seconds += int(match.group(1)) * 3600
    match = re.search(r'(\d+)m(?!s)', duration_str)
    if match:
        total_seconds += int(match.group(1)) * 60
    match = re.search(r'([\d.]+)s', duration_str)
    if match:
        total_seconds += float(match.group(1))
    return total_seconds if total_seconds > 0 else 0


def parse_memory(mem_str):
    """Parse memory string to bytes."""
    if not mem_str or mem_str == '-':
        return 0
    mem_str = str(mem_str).strip().upper()
    try:
        return int(float(mem_str))
    except ValueError:
        pass
    units = {'B': 1, 'KB': 1024, 'MB': 1024**2, 'GB': 1024**3, 'TB': 1024**4,
             'K': 1024, 'M': 1024**2, 'G': 1024**3, 'T': 1024**4}
    match = re.match(r'([\d.]+)\s*([KMGT]?B?)', mem_str)
    if match:
        value = float(match.group(1))
        unit = match.group(2) or 'B'
        return int(value * units.get(unit, 1))
    return 0


def format_bytes(bytes_val):
    """Format bytes to human-readable string."""
    if bytes_val == 0:
        return "0 B"
    units = ['B', 'KB', 'MB', 'GB', 'TB']
    i = 0
    val = float(bytes_val)
    while val >= 1024 and i < len(units) - 1:
        val /= 1024
        i += 1
    return f"{val:.1f} {units[i]}"


def format_duration_str(seconds):
    """Format seconds to human-readable duration string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = seconds % 60
        return f"{mins}m {secs:.0f}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


def parse_timestamp(ts_str):
    """Parse Nextflow timestamp string to datetime."""
    if not ts_str or ts_str == '-':
        return None
    for fmt in ['%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S']:
        try:
            return datetime.strptime(ts_str, fmt)
        except ValueError:
            continue
    return None


def parse_newick(newick_str):
    """Parse a Newick format string into a tree structure.

    Returns a dict with 'name', 'children', and 'length' keys.
    """
    idx = 0

    def parse_node():
        nonlocal idx
        node = {'name': '', 'children': [], 'length': 0}

        if newick_str[idx] == '(':
            idx += 1  # skip '('
            while newick_str[idx] != ')':
                node['children'].append(parse_node())
                if newick_str[idx] == ',':
                    idx += 1
            idx += 1  # skip ')'

        # Parse name and branch length
        name_end = idx
        while name_end < len(newick_str) and newick_str[name_end] not in ',):;':
            name_end += 1

        name_part = newick_str[idx:name_end]
        if ':' in name_part:
            name, length = name_part.rsplit(':', 1)
            node['name'] = name.strip("'\"")
            try:
                node['length'] = float(length)
            except ValueError:
                node['length'] = 0
        else:
            node['name'] = name_part.strip("'\"")

        idx = name_end
        return node

    return parse_node()


def sanitize_taxon(name):
    """Sanitize taxon name for use in file paths."""
    result = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    result = re.sub(r'_+', '_', result)
    return result.strip('_')
