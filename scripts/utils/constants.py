#!/usr/bin/env python3
"""Shared constants for BGC analysis pipeline."""

# BGC type colors (matches antiSMASH conventions)
BGC_COLORS = {
    'NRPS': '#1f77b4',
    'T1PKS': '#ff7f0e',
    'T2PKS': '#2ca02c',
    'T3PKS': '#d62728',
    'terpene': '#9467bd',
    'RiPP': '#8c564b',
    'lanthipeptide': '#e377c2',
    'bacteriocin': '#7f7f7f',
    'siderophore': '#bcbd22',
    'phosphonate': '#17becf',
    'other': '#aec7e8',
    'hybrid': '#ff9896',
    'NAPAA': '#98df8a',
    'NRP-metallophore': '#c5b0d5',
    'ladderane': '#c49c94',
    'butyrolactone': '#f7b6d2',
    'ectoine': '#c7c7c7',
    'NAGGN': '#dbdb8d',
    'hserlactone': '#9edae5',
    'arylpolyene': '#393b79',
    'resorcinol': '#637939',
    'phenazine': '#8c6d31',
    'melanin': '#843c39',
    'betalactone': '#7b4173',
}

# Gene function colors (matches antiSMASH conventions)
GENE_COLORS = {
    'biosynthetic': '#e74c3c',           # Red - core biosynthetic
    'biosynthetic-additional': '#e67e22', # Orange - additional biosynthetic
    'regulatory': '#27ae60',              # Green - regulatory
    'transport': '#3498db',               # Blue - transport
    'resistance': '#9b59b6',              # Purple - resistance
    'other': '#95a5a6',                   # Gray - other/unknown
}

# GCF color palette for tree visualization
GCF_COLORS = [
    '#e41a1c',  # Red
    '#377eb8',  # Blue
    '#4daf4a',  # Green
    '#984ea3',  # Purple
    '#ff7f00',  # Orange
    '#ffff33',  # Yellow
    '#a65628',  # Brown
    '#f781bf',  # Pink
    '#999999',  # Gray
    '#66c2a5',  # Teal
    '#fc8d62',  # Coral
    '#8da0cb',  # Light blue
    '#e78ac3',  # Magenta
    '#a6d854',  # Lime
    '#ffd92f',  # Gold
    '#e5c494',  # Tan
]

# KCB similarity thresholds
KCB_THRESHOLDS = {
    'high': 75,
    'medium': 50,
    'low': 15,
}

# KCB similarity colors
KCB_COLORS = {
    'high': '#28a745',    # Green
    'medium': '#fd7e14',  # Orange
    'low': '#6c757d',     # Gray
}
