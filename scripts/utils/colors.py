"""Shared color utility functions."""

import colorsys
import hashlib


def rank_or_hash_color(seed_str, rank, n_top,
                        vivid_sat=0.75, vivid_val=0.85,
                        muted_sat=0.35, muted_val=0.75):
    """
    Assign a hex color based on rank within n_top.
    - rank < n_top : evenly-spaced vivid hue
    - rank >= n_top: muted hash-derived hue (deterministic per seed_str)
    """
    if rank < n_top:
        hue = rank / n_top
        r, g, b = colorsys.hsv_to_rgb(hue, vivid_sat, vivid_val)
    else:
        h = int(hashlib.md5(seed_str.encode()).hexdigest(), 16)
        hue = (h % 360) / 360.0
        r, g, b = colorsys.hsv_to_rgb(hue, muted_sat, muted_val)
    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))


def family_color(family_id, rank, n_top):
    """Assign a color to a GCF family by size rank. None → light gray."""
    if family_id is None:
        return '#cccccc'
    return rank_or_hash_color(str(family_id), rank, n_top, vivid_sat=0.75)


def genus_color(genus, rank, n_top):
    """Assign a color to a genus by frequency rank."""
    return rank_or_hash_color(genus, rank, n_top, vivid_sat=0.7)
