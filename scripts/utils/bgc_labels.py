"""Shared BGC label utilities."""

import os


def label_from_path(gbk_path):
    """
    Derive a short BGC label from a GBK file path.
    Path pattern: .../antismash_results/{taxon}/{genome}/{genome}.region{NNN}.gbk
    Returns '{genome}.region{NNN}' truncated to 60 chars if needed.
    Sanitizes characters that confuse Newick parsers (commas, colons, parens).
    """
    basename = os.path.basename(gbk_path)
    name = basename.replace('.gbk', '').replace('.gbff', '')
    name = name.replace(',', '_').replace(':', '_').replace('(', '_').replace(')', '_')
    return name[:60] if len(name) > 60 else name


def make_labels_unique(labels):
    """Append a numeric suffix to any duplicate labels."""
    seen = {}
    result = []
    for label in labels:
        if label not in seen:
            seen[label] = 0
            result.append(label)
        else:
            seen[label] += 1
            result.append(f'{label}_{seen[label]}')
    return result
