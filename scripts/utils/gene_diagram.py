#!/usr/bin/env python3
"""Gene diagram SVG generation for BGC visualization."""


def generate_gene_svg(genes, region_start, region_end, width=900, height=80):
    """Generate SVG gene diagram for a BGC with proper arrow shapes.

    Args:
        genes: List of gene dicts with keys: start, end, strand, color, locus_tag, product, gene_name
        region_start: Start position of the region
        region_end: End position of the region
        width: SVG width in pixels
        height: SVG height in pixels

    Returns:
        SVG string for the gene diagram
    """
    if not genes:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="80"><text x="450" y="40" text-anchor="middle" fill="#999">No gene data available</text></svg>'

    region_length = region_end - region_start
    if region_length <= 0:
        region_length = 1

    scale = (width - 100) / region_length  # Leave margins
    margin_left = 50
    gene_y = height / 2
    arrow_height = 20
    arrow_tip_width = 12  # Width of the arrowhead

    svg_parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" style="background: #f8f9fa;">',
        # Scale bar
        f'<line x1="{margin_left}" y1="{height - 10}" x2="{width - margin_left}" y2="{height - 10}" stroke="#ccc" stroke-width="1"/>',
        f'<text x="{margin_left}" y="{height - 2}" font-size="9" fill="#666">{region_start:,}</text>',
        f'<text x="{width - margin_left}" y="{height - 2}" font-size="9" fill="#666" text-anchor="end">{region_end:,}</text>',
    ]

    for gene in genes:
        gene_start = (gene['start'] - region_start) * scale + margin_left
        gene_end = (gene['end'] - region_start) * scale + margin_left
        gene_width = max(gene_end - gene_start, 8)  # Minimum width

        color = gene['color']
        half_height = arrow_height / 2
        body_half_height = arrow_height / 3  # Body is narrower than arrowhead

        # Tooltip text
        tooltip = f"{gene['locus_tag']}: {gene['product']}"
        if gene.get('gene_name'):
            tooltip = f"{gene['gene_name']} ({gene['locus_tag']}): {gene['product']}"

        if gene['strand'] == 1:  # Forward strand (right-pointing arrow)
            if gene_width > arrow_tip_width:
                # Arrow with body: 7-point polygon
                # Shape: rectangle body + triangular arrowhead pointing right
                body_end = gene_end - arrow_tip_width
                points = (
                    f"{gene_start},{gene_y - body_half_height} "  # Top-left of body
                    f"{body_end},{gene_y - body_half_height} "    # Top-right of body
                    f"{body_end},{gene_y - half_height} "         # Top of arrowhead base
                    f"{gene_end},{gene_y} "                       # Arrow tip
                    f"{body_end},{gene_y + half_height} "         # Bottom of arrowhead base
                    f"{body_end},{gene_y + body_half_height} "    # Bottom-right of body
                    f"{gene_start},{gene_y + body_half_height}"   # Bottom-left of body
                )
            else:
                # Small gene: simple triangle pointing right
                points = (
                    f"{gene_start},{gene_y - half_height} "
                    f"{gene_end},{gene_y} "
                    f"{gene_start},{gene_y + half_height}"
                )
        else:  # Reverse strand (left-pointing arrow)
            if gene_width > arrow_tip_width:
                # Arrow with body: 7-point polygon
                # Shape: triangular arrowhead pointing left + rectangle body
                body_start = gene_start + arrow_tip_width
                points = (
                    f"{gene_start},{gene_y} "                     # Arrow tip
                    f"{body_start},{gene_y - half_height} "       # Top of arrowhead base
                    f"{body_start},{gene_y - body_half_height} "  # Top-left of body
                    f"{gene_end},{gene_y - body_half_height} "    # Top-right of body
                    f"{gene_end},{gene_y + body_half_height} "    # Bottom-right of body
                    f"{body_start},{gene_y + body_half_height} "  # Bottom-left of body
                    f"{body_start},{gene_y + half_height}"        # Bottom of arrowhead base
                )
            else:
                # Small gene: simple triangle pointing left
                points = (
                    f"{gene_end},{gene_y - half_height} "
                    f"{gene_start},{gene_y} "
                    f"{gene_end},{gene_y + half_height}"
                )

        svg_parts.append(
            f'<polygon points="{points}" fill="{color}" stroke="#333" stroke-width="0.5" opacity="0.9">'
            f'<title>{tooltip}</title>'
            f'</polygon>'
        )

    svg_parts.append('</svg>')
    return ''.join(svg_parts)
