#!/usr/bin/env python3
"""Clustering visualization functions for BiG-SCAPE and BiG-SLiCE results."""

import json
import os
import sys
import urllib.parse
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def generate_bigslice_stats_html(bigslice_stats_file):
    """Generate HTML for BiG-SLiCE clustering statistics."""
    if not os.path.exists(bigslice_stats_file) or os.path.basename(bigslice_stats_file).startswith('NO_'):
        return ''

    try:
        with open(bigslice_stats_file, 'r') as f:
            stats = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read BiG-SLiCE statistics: {e}")
        return ''

    if 'error' in stats:
        return f'''
    <h3>BiG-SLiCE Clustering Statistics</h3>
    <div class="info-box" style="background-color: #fff3cd; border-left: 4px solid #ffc107;">
        <p><strong>Note:</strong> {stats['error']}</p>
    </div>
    '''

    # Generate fragmentation breakdown
    frag_html = ''
    if 'fragmentation' in stats:
        frag_rows = []
        for frag_type, frag_data in stats['fragmentation'].items():
            frag_rows.append(f'''
                <tr>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd;">{frag_type.capitalize()}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{frag_data['count']}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{frag_data['avg_distance']}</td>
                </tr>
            ''')
        frag_html = f'''
            <tr>
                <td colspan="3" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">Fragmentation Breakdown:</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold;">Type</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">Count</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">Avg Distance</td>
            </tr>
            {''.join(frag_rows)}
        '''

    # Generate GCF size distribution preview (top 5)
    gcf_dist_html = ''
    if 'gcf_distribution' in stats and stats['gcf_distribution']:
        top_gcfs = stats['gcf_distribution'][:5]
        gcf_rows = []
        for gcf in top_gcfs:
            gcf_rows.append(f'''
                <tr>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd;">GCF {gcf['gcf_id']}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{gcf['bgc_count']}</td>
                </tr>
            ''')
        more_text = ''
        if len(stats['gcf_distribution']) > 5:
            more_text = f'''
                <tr>
                    <td colspan="2" style="padding: 8px; text-align: center; color: #666; font-style: italic;">
                        ... and {len(stats['gcf_distribution']) - 5} more GCFs
                    </td>
                </tr>
            '''
        gcf_dist_html = f'''
            <tr>
                <td colspan="2" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">Top Gene Cluster Families:</td>
            </tr>
            {''.join(gcf_rows)}
            {more_text}
        '''

    return f'''
    <h3>BiG-SLiCE Clustering Statistics</h3>
    <div class="info-box">
        <table style="width: 100%; border-collapse: collapse;">
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; font-weight: bold;">Total BGCs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_bgcs', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; font-weight: bold;">Total GCFs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_gcfs', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Distance Threshold (T)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('threshold', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Avg. Distance to GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('avg_distance_to_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">BGCs Assigned (d <= T)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('bgcs_assigned', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">BGCs Not Assigned (d > T)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('bgcs_not_assigned', 0)}</td>
            </tr>
            <tr>
                <td colspan="2" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">GCF Size Statistics:</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Min BGCs per GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('min_bgcs_per_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Avg BGCs per GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('avg_bgcs_per_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Max BGCs per GCF</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('max_bgcs_per_gcf', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Singleton GCFs (n=1)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('singleton_gcfs', 0)}</td>
            </tr>
            {frag_html}
            {gcf_dist_html}
        </table>
        <p style="margin-top: 15px; font-size: 0.9em; color: #666;">
            <strong>Note:</strong> Statistics extracted from BiG-SLiCE database. For detailed analysis, use the BiG-SLiCE web interface or query the database directly.
        </p>
    </div>
    '''


def generate_bigscape_stats_html(bigscape_stats_file, mibig_included=False):
    """Generate HTML for BiG-SCAPE clustering statistics."""
    if not os.path.exists(bigscape_stats_file) or os.path.basename(bigscape_stats_file).startswith('NO_'):
        return ''

    try:
        with open(bigscape_stats_file, 'r') as f:
            stats = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read BiG-SCAPE statistics: {e}")
        return ''

    if 'error' in stats:
        return f'''
    <h3>BiG-SCAPE Clustering Statistics</h3>
    <div class="info-box" style="background-color: #fff3cd; border-left: 4px solid #ffc107;">
        <p><strong>Note:</strong> {stats['error']}</p>
    </div>
    '''

    # Generate BGC class breakdown if available
    class_breakdown_html = ''
    if 'bgc_classes' in stats and stats['bgc_classes']:
        class_rows = []
        for bgc_class, class_data in sorted(stats['bgc_classes'].items()):
            class_rows.append(f'''
                <tr>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd;">{bgc_class}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{class_data['families']}</td>
                    <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{class_data['bgcs']}</td>
                </tr>
            ''')
        class_breakdown_html = f'''
            <tr>
                <td colspan="2" style="padding: 12px 8px 4px 8px; font-weight: bold; color: #666;">BGC Class Breakdown:</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold;">Class</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">Families</td>
                <td style="padding: 8px; border-bottom: 2px solid #999; font-weight: bold; text-align: right;">BGCs</td>
            </tr>
            {''.join(class_rows)}
        '''

    # Build the complete statistics table
    cutoff_display = stats.get('cutoff', 'N/A')

    # Show N/A for MIBiG families if MIBiG was not included in the analysis
    if mibig_included:
        mibig_families_display = stats.get('families_with_mibig', 0)
    else:
        mibig_families_display = '<span style="color: #999;" title="MIBiG references were not included in this analysis">N/A</span>'

    return f'''
    <h3>BiG-SCAPE Clustering Statistics</h3>
    <div class="info-box">
        <table style="width: 100%; border-collapse: collapse;">
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; font-weight: bold;">Network Overview</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;"></td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Number of Families</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_families', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Total BGCs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right; font-weight: bold;">{stats.get('total_bgcs', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Average BGCs per Family</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('avg_bgcs_per_family', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Max BGCs in a Family</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('max_bgcs_per_family', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Singleton Families (n=1)</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{stats.get('singleton_families', 0)}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Families with MIBiG Reference BGCs</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{mibig_families_display}</td>
            </tr>
            <tr>
                <td style="padding: 8px; border-bottom: 1px solid #ddd;">Distance Cutoff</td>
                <td style="padding: 8px; border-bottom: 1px solid #ddd; text-align: right;">{cutoff_display}</td>
            </tr>
            {class_breakdown_html}
        </table>
        <p style="margin-top: 15px; font-size: 0.9em; color: #666;">
            <strong>Note:</strong> Statistics extracted from BiG-SCAPE clustering files. For detailed network analysis, open the BiG-SCAPE HTML output in <code>results/bigscape_results/[taxon]/index.html</code>.
        </p>
    </div>
    '''


def generate_gcf_visualization_html(gcf_data_file, taxon):
    """Generate HTML for GCF representative visualization in Clustering tab."""
    if not os.path.exists(gcf_data_file):
        return ''

    try:
        with open(gcf_data_file, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Warning: Could not read GCF data: {e}")
        return ''

    if 'error' in data:
        return f'''
    <div class="gcf-section">
        <h4>Gene Cluster Family Representatives</h4>
        <div class="info-box" style="background-color: #fff3cd; border-left: 4px solid #ffc107;">
            <p><strong>Note:</strong> {data['error']}</p>
        </div>
    </div>
    '''

    gcfs = data.get('gcfs', [])
    summary = data.get('summary', {})

    if not gcfs:
        return '''
    <div class="gcf-section">
        <h4>Gene Cluster Family Representatives</h4>
        <p style="color: #666;">No gene cluster families found.</p>
    </div>
    '''

    total = summary.get('total', len(gcfs))
    singletons = summary.get('singletons', 0)
    clusters = summary.get('clusters', 0)

    # Build GCF cards HTML
    gcf_cards = ''
    for gcf in gcfs:
        family_id = gcf.get('family_id', '?')
        member_count = gcf.get('member_count', 0)
        is_singleton = gcf.get('is_singleton', False)
        product = gcf.get('product', 'Unknown')
        organism = gcf.get('organism', 'Unknown organism')
        genome_name = gcf.get('genome_name', 'unknown')
        region_number = gcf.get('region_number', 1)
        svg_diagram = gcf.get('svg_diagram', '')
        genes = gcf.get('genes', [])
        antismash_link = gcf.get('antismash_link', '')

        badge_class = 'singleton' if is_singleton else 'cluster'
        badge_text = 'Singleton' if is_singleton else f'{member_count} BGCs'

        # Check if domain hits are stored directly in genes (new approach)
        # or in enhanced_analysis (legacy approach)
        has_enhanced = False
        use_gene_level_hits = False

        # Check first gene for direct domain hits
        if genes and (genes[0].get('clusterhmmer_hits') is not None or genes[0].get('tigrfam_hits') is not None):
            use_gene_level_hits = True
            has_enhanced = any(
                gene.get('clusterhmmer_hits') or gene.get('tigrfam_hits')
                for gene in genes
            )

        # Fallback: build lookup dicts from enhanced_analysis (legacy)
        clusterhmmer_by_gene = {}
        tigrfam_by_gene = {}
        if not use_gene_level_hits:
            enhanced_analysis = gcf.get('enhanced_analysis')
            if enhanced_analysis:
                clusterhmmer_hits = enhanced_analysis.get('clusterhmmer_hits', [])
                tigrfam_hits = enhanced_analysis.get('tigrfam_hits', [])
                has_enhanced = bool(clusterhmmer_hits or tigrfam_hits)

                # Build lookup by gene for ClusterHmmer
                for hit in clusterhmmer_hits:
                    gene_id = hit.get('gene', '')
                    if gene_id:
                        if gene_id not in clusterhmmer_by_gene:
                            clusterhmmer_by_gene[gene_id] = []
                        clusterhmmer_by_gene[gene_id].append(hit)

                # Build lookup by gene for TIGRFam
                for hit in tigrfam_hits:
                    gene_id = hit.get('gene', '')
                    if gene_id:
                        if gene_id not in tigrfam_by_gene:
                            tigrfam_by_gene[gene_id] = []
                        tigrfam_by_gene[gene_id].append(hit)

        # Build gene table rows with integrated enhanced analysis
        gene_rows = ''
        for gene in genes[:30]:  # Limit to first 30 genes for display
            locus_tag = gene.get('locus_tag', '-')
            gene_name = gene.get('gene_name', '-') or '-'
            gene_product = gene.get('product', 'hypothetical protein')
            function = gene.get('function', 'other')
            color = gene.get('color', '#95a5a6')

            # Get ClusterHmmer domains for this gene
            if use_gene_level_hits:
                ch_hits = gene.get('clusterhmmer_hits', [])
            else:
                ch_hits = clusterhmmer_by_gene.get(locus_tag, [])
            # Get TIGRFam annotations for this gene
            if use_gene_level_hits:
                tf_hits = gene.get('tigrfam_hits', [])
            else:
                tf_hits = tigrfam_by_gene.get(locus_tag, [])

            # Enhance product annotation if it's generic
            # Priority: TIGRFam description > Pfam description > original
            generic_products = ['hypothetical protein', 'putative protein', 'unknown',
                               'uncharacterized protein', 'predicted protein', 'conserved protein',
                               'domain-containing protein', 'family protein']
            product_is_generic = any(gp in gene_product.lower() for gp in generic_products)
            enhanced_product = gene_product
            product_source = ''

            if product_is_generic or gene_product == '-':
                # Try TIGRFam first (more curated)
                if tf_hits:
                    best_tf = max(tf_hits, key=lambda x: x.get('score', 0))
                    tf_desc = best_tf.get('description', '')
                    if tf_desc and len(tf_desc) > 5:
                        # Clean up TIGRFam description (remove prefix like "TIGR02320: ")
                        if ': ' in tf_desc:
                            tf_desc = tf_desc.split(': ', 1)[-1]
                        enhanced_product = tf_desc
                        product_source = 'TF'
                # Fall back to Pfam if no TIGRFam
                if product_source == '' and ch_hits:
                    best_ch = max(ch_hits, key=lambda x: x.get('score', 0))
                    ch_desc = best_ch.get('description', '')
                    if ch_desc and len(ch_desc) > 5:
                        enhanced_product = ch_desc
                        product_source = 'Pfam'

            # Format product cell with source indicator
            if product_source:
                product_cell = f'{enhanced_product[:35]}{"..." if len(enhanced_product) > 35 else ""} <span style="background: #3498db; color: white; padding: 1px 4px; border-radius: 3px; font-size: 0.7em;" title="Inferred from {product_source}">{product_source}</span>'
            else:
                product_cell = f'{gene_product[:40]}{"..." if len(gene_product) > 40 else ""}'

            # Format ClusterHmmer column
            if ch_hits:
                ch_domains = [f"{h.get('hmm', '')} ({h.get('score', 0):.0f})" for h in ch_hits[:3]]
                clusterhmmer_col = '<br>'.join(ch_domains)
                if len(ch_hits) > 3:
                    clusterhmmer_col += f'<br><span style="color: #888;">+{len(ch_hits)-3} more</span>'
            else:
                clusterhmmer_col = '<span style="color: #ccc;">-</span>'

            # Format TIGRFam column
            if tf_hits:
                tf_annots = [f"{h.get('domain', '')} ({h.get('score', 0):.0f})" for h in tf_hits[:2]]
                tigrfam_col = '<br>'.join(tf_annots)
                if len(tf_hits) > 2:
                    tigrfam_col += f'<br><span style="color: #888;">+{len(tf_hits)-2} more</span>'
            else:
                tigrfam_col = '<span style="color: #ccc;">-</span>'

            # Format BlastP link
            translation = gene.get('translation', '')
            if translation:
                blastp_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&QUERY={urllib.parse.quote(translation)}&LINK_LOC=protein&PAGE_TYPE=BlastSearch"
                blastp_col = f'<a href="{blastp_url}" target="_blank" style="display: inline-block; background: #2c5aa0; color: white; padding: 2px 6px; border-radius: 3px; font-size: 0.75em; text-decoration: none;" title="Search NCBI BlastP">BlastP</a>'
            else:
                blastp_col = '<span style="color: #ccc;">-</span>'

            gene_rows += f'''
                <tr>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">
                        <span style="display: inline-block; width: 10px; height: 10px; background: {color}; border-radius: 2px; margin-right: 6px;"></span>
                        {locus_tag}
                    </td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">{gene_name}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">{product_cell}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee;">{function}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee; font-size: 0.8em;">{clusterhmmer_col}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee; font-size: 0.8em;">{tigrfam_col}</td>
                    <td style="padding: 6px 8px; border-bottom: 1px solid #eee; text-align: center;">{blastp_col}</td>
                </tr>
            '''

        more_genes_note = f'<p style="color: #666; font-size: 0.85em; margin-top: 5px;">Showing first 30 of {len(genes)} genes</p>' if len(genes) > 30 else ''

        # AntiSMASH link
        antismash_link_html = f'<a href="{antismash_link}" target="_blank" style="color: #2c5aa0;">View in antiSMASH</a>' if antismash_link else ''

        gcf_cards += f'''
        <div class="gcf-card" data-type="{badge_class}" data-size="{member_count}" data-product="{product}">
            <div class="gcf-header" onclick="toggleGCF('gcf_{family_id}')">
                <span class="gcf-title">GCF {family_id}: {product}</span>
                <span class="gcf-badge {badge_class}">{badge_text}</span>
                <span class="gcf-organism">{organism}</span>
                <span class="gcf-toggle" id="gcf_{family_id}_toggle">+</span>
            </div>
            <div class="gcf-content" id="gcf_{family_id}_content" style="display: none;">
                <div class="gene-diagram">
                    {svg_diagram if svg_diagram else '<p style="color: #999;">Gene diagram not available</p>'}
                </div>
                <div class="gene-legend">
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #e74c3c; border-radius: 2px;"></span> Core biosynthetic</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #e67e22; border-radius: 2px;"></span> Additional biosynthetic</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #27ae60; border-radius: 2px;"></span> Regulatory</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #3498db; border-radius: 2px;"></span> Transport</span>
                    <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #9b59b6; border-radius: 2px;"></span> Resistance</span>
                    <span><span style="display: inline-block; width: 12px; height: 12px; background: #95a5a6; border-radius: 2px;"></span> Other</span>
                </div>
                <div class="gene-table-container" style="margin-top: 15px; max-height: 400px; overflow-y: auto;">
                    <table class="gene-table" style="width: 100%; border-collapse: collapse; font-size: 0.85em;">
                        <thead>
                            <tr style="background: #f5f5f5;">
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Locus Tag</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Gene</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Product</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Function</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">ClusterHmmer</th>
                                <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">TIGRFam</th>
                                <th style="padding: 8px; text-align: center; border-bottom: 2px solid #ddd;">BlastP</th>
                            </tr>
                        </thead>
                        <tbody>
                            {gene_rows}
                        </tbody>
                    </table>
                    {more_genes_note}
                </div>
                <div class="gcf-links" style="margin-top: 15px; padding-top: 10px; border-top: 1px solid #eee;">
                    {antismash_link_html}
                </div>
            </div>
        </div>
        '''

    # Build complete GCF visualization HTML with styles and scripts
    gcf_html = f'''
    <div class="gcf-section" style="margin-top: 30px;">
        <h4>Gene Cluster Family Representatives</h4>
        <p style="color: #666; margin-bottom: 15px;">
            Representative BGCs for each Gene Cluster Family identified by BiG-SCAPE. Click on a family to expand and view the gene diagram and annotations.
        </p>

        <div class="gcf-summary" style="display: flex; gap: 15px; margin-bottom: 20px;">
            <div style="background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center;">
                <div style="font-size: 1.5em; font-weight: bold; color: #2c5aa0;">{total}</div>
                <div style="font-size: 0.85em; color: #666;">Total GCFs</div>
            </div>
            <div style="background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center;">
                <div style="font-size: 1.5em; font-weight: bold; color: #27ae60;">{clusters}</div>
                <div style="font-size: 0.85em; color: #666;">Clusters (2+ BGCs)</div>
            </div>
            <div style="background: #f8f9fa; padding: 12px 20px; border-radius: 8px; text-align: center;">
                <div style="font-size: 1.5em; font-weight: bold; color: #f39c12;">{singletons}</div>
                <div style="font-size: 0.85em; color: #666;">Singletons</div>
            </div>
        </div>

        <div class="gcf-controls" style="margin-bottom: 15px; display: flex; gap: 10px; align-items: center;">
            <label style="font-size: 0.9em; color: #666;">Filter:</label>
            <select id="gcfFilter" onchange="filterGCFs(this.value)" style="padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px;">
                <option value="all">All GCFs</option>
                <option value="cluster">Clusters Only</option>
                <option value="singleton">Singletons Only</option>
            </select>
            <label style="font-size: 0.9em; color: #666; margin-left: 15px;">Sort:</label>
            <select id="gcfSort" onchange="sortGCFs(this.value)" style="padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px;">
                <option value="size">By Size (largest first)</option>
                <option value="product">By Product Type</option>
            </select>
            <button onclick="toggleAllGCFs()" style="margin-left: 15px; padding: 6px 12px; border: 1px solid #ddd; border-radius: 4px; background: #f8f9fa; cursor: pointer;">Expand/Collapse All</button>
        </div>

        <div class="gcf-container" id="gcfContainer">
            {gcf_cards}
        </div>
    </div>

    <style>
        .gcf-card {{
            border: 1px solid #ddd;
            border-radius: 8px;
            margin-bottom: 10px;
            background: white;
            overflow: hidden;
        }}
        .gcf-header {{
            padding: 12px 15px;
            cursor: pointer;
            display: flex;
            align-items: center;
            gap: 12px;
            background: #fafafa;
            transition: background 0.2s;
        }}
        .gcf-header:hover {{
            background: #f0f0f0;
        }}
        .gcf-title {{
            font-weight: 600;
            color: #333;
            flex-grow: 1;
        }}
        .gcf-badge {{
            padding: 3px 10px;
            border-radius: 12px;
            font-size: 0.8em;
            font-weight: 500;
        }}
        .gcf-badge.singleton {{
            background: #fff3cd;
            color: #856404;
        }}
        .gcf-badge.cluster {{
            background: #d4edda;
            color: #155724;
        }}
        .gcf-organism {{
            color: #666;
            font-size: 0.85em;
            font-style: italic;
        }}
        .gcf-toggle {{
            font-size: 1.2em;
            color: #999;
            font-weight: bold;
        }}
        .gcf-content {{
            padding: 15px;
            border-top: 1px solid #eee;
        }}
        .gene-diagram {{
            padding: 10px;
            background: #f8f9fa;
            border-radius: 4px;
            overflow-x: auto;
        }}
        .gene-legend {{
            margin-top: 10px;
            font-size: 0.8em;
            color: #666;
        }}
    </style>

    <script>
        function toggleGCF(gcfId) {{
            const content = document.getElementById(gcfId + '_content');
            const toggle = document.getElementById(gcfId + '_toggle');
            if (content.style.display === 'none') {{
                content.style.display = 'block';
                toggle.textContent = '-';
            }} else {{
                content.style.display = 'none';
                toggle.textContent = '+';
            }}
        }}

        function filterGCFs(type) {{
            const cards = document.querySelectorAll('.gcf-card');
            cards.forEach(card => {{
                if (type === 'all' || card.dataset.type === type) {{
                    card.style.display = 'block';
                }} else {{
                    card.style.display = 'none';
                }}
            }});
        }}

        function sortGCFs(by) {{
            const container = document.getElementById('gcfContainer');
            const cards = Array.from(container.querySelectorAll('.gcf-card'));
            cards.sort((a, b) => {{
                if (by === 'size') {{
                    return parseInt(b.dataset.size) - parseInt(a.dataset.size);
                }}
                return a.dataset.product.localeCompare(b.dataset.product);
            }});
            cards.forEach(card => container.appendChild(card));
        }}

        function toggleAllGCFs() {{
            const contents = document.querySelectorAll('.gcf-content');
            const anyHidden = Array.from(contents).some(c => c.style.display === 'none');
            contents.forEach(content => {{
                content.style.display = anyHidden ? 'block' : 'none';
                const id = content.id.replace('_content', '_toggle');
                const toggle = document.getElementById(id);
                if (toggle) toggle.textContent = anyHidden ? '-' : '+';
            }});
        }}
    </script>
    '''

    return gcf_html
