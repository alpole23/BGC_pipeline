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

# Coupling enzyme class colors (used across bgc_gcf_tree, bgc_gcf_heatmap,
# bgc_all_bgcs_tree, bgc_coupling_tree, bgc_coupling_annotation)
COUPLING_COLORS = {
    'FrbC-like': '#e41a1c',
    'Ppd':       '#377eb8',
    'Ppd-CDP':   '#984ea3',
    'VlpB-like': '#4daf4a',
    'PalB-like': '#ff7f00',
    'Unknown':   '#aaaaaa',
}

# Display order for coupling enzyme classes in legends
COUPLING_ORDER = ['FrbC-like', 'VlpB-like', 'Ppd', 'Ppd-CDP', 'PalB-like', 'Unknown']

# Normalize legacy coupling class names from older annotation files
LEGACY_CLASS_NAMES = {
    'Fe-ADH':  'VlpB-like',
    'TPP+NTP': 'Ppd-CDP',
    'PalB':    'PalB-like',
    'FrbC':    'FrbC-like',
}

# Pfam accession → short human-readable name.
# Verified against antiSMASH clusterhmmer output on phosphonate BGCs.
# Add entries here as needed; unknown domains fall back to their accession.
DOMAIN_NAMES = {
    'PF13714': 'PEP_mutase',
    'PF00296': 'HEPD',           # 2-hydroxyethylphosphonate dioxygenase
    'PF02775': 'ThDP_C',         # phosphonopyruvate decarboxylase (C-term)
    'PF02776': 'ThDP_N',         # phosphonopyruvate decarboxylase (N-term)
    'PF00266': 'Aminotrans_V',   # 2-AEP transaminase (class V)
    'PF00155': 'Aminotrans_I',   # aminotransferase class I/II
    'PF00682': 'FrbC_HMGL',      # FrbC-like / phosphonomethylmalate synthase (HMGL superfamily)
    'PF13649': 'Radical_SAM',
    'PF08241': 'Methyltransf',
    'PF00330': 'Aconitase',
    'PF00694': 'Aconitase_C',
    'PF00149': 'Metallophos',    # calcineurin-like phosphoesterase
    'PF00571': 'CBS',
    'PF07690': 'MFS_transporter',
    'PF00005': 'ABC_ATPase',
    'PF01118': 'SemiAldDH_NAD',  # semialdehyde dehydrogenase
    'PF02774': 'SemiAldDH_C',
    'PF00464': 'SHMT',           # serine hydroxymethyltransferase
    'PF12804': 'NTP_transf',
    'PF13673': 'Acetyltransf',
    'PF07228': 'DUF1453',
    'PF00733': 'Asn_synthase',
    'PF13537': 'GATase',
    'PF00581': 'Rhodanese',
    'PF01613': 'Trp_syntA',
    'PF00202': 'Aminotrans_III',
    'PF00892': 'EamA_transporter',
    'PF22617': 'PF22617',
    'PF05321': 'PF05321',
}
