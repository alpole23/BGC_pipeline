#!/usr/bin/env python3
"""Taxonomy tree visualization functions."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def generate_taxonomy_tree_html(taxonomy_tree_data):
    """Generate interactive HTML for the taxonomy tree with genome-level data."""
    tree = taxonomy_tree_data.get('tree', {})
    metadata = taxonomy_tree_data.get('metadata', {})

    def render_genome_list(genomes, level):
        """Render genome list under a species node with heatmap coloring."""
        if not genomes:
            return ''

        # Get all unique BGC types from genomes
        bgc_type_columns = set()
        for genome_data in genomes:
            bgc_types = genome_data.get('bgc_types', {})
            bgc_type_columns.update(bgc_types.keys())
        bgc_type_columns = sorted(bgc_type_columns)

        # Calculate global max count across ALL BGC types for consistent heatmap coloring
        global_bgc_max = 0
        for genome_data in genomes:
            bgc_types = genome_data.get('bgc_types', {})
            for count in bgc_types.values():
                if count > global_bgc_max:
                    global_bgc_max = count
        global_bgc_max = global_bgc_max if global_bgc_max > 0 else 1

        # Calculate max total BGCs for total column heatmap
        total_max = max([g.get('total_bgcs', 0) for g in genomes])
        total_max = total_max if total_max > 0 else 1

        def get_green_bg_color(count, max_count):
            """Convert count to sequential green background color for Total BGCs."""
            if count == 0:
                return ''
            intensity = count / max_count
            if intensity <= 0.2:
                return '#e8f5e9'
            elif intensity <= 0.4:
                return '#a5d6a7'
            elif intensity <= 0.6:
                return '#66bb6a'
            elif intensity <= 0.8:
                return '#43a047'
            else:
                return '#2e7d32'

        def get_red_font_color(count, max_count):
            """Convert count to sequential red font color for BGC types."""
            if count == 0:
                return ''
            intensity = count / max_count
            if intensity <= 0.2:
                return '#ffcdd2'
            elif intensity <= 0.4:
                return '#ef5350'
            elif intensity <= 0.6:
                return '#e53935'
            elif intensity <= 0.8:
                return '#c62828'
            else:
                return '#b71c1c'

        genome_list_html = '<div class="genome-list">'
        genome_list_html += '<table class="genome-table">'
        genome_list_html += '<tr><th>Genome</th><th>Total BGCs</th>'

        for bgc_type in bgc_type_columns:
            genome_list_html += f'<th>{bgc_type}</th>'
        genome_list_html += '</tr>'

        sorted_genomes = sorted(genomes, key=lambda x: x.get('name', ''))

        for genome_data in sorted_genomes:
            genome_name = genome_data.get('name', 'Unknown')
            total = genome_data.get('total_bgcs', 0)
            bgc_types = genome_data.get('bgc_types', {})

            genome_list_html += f'<tr><td><a href="genomes/{genome_name}.html" class="genome-link">{genome_name}</a></td>'

            total_bg_color = get_green_bg_color(total, total_max)
            total_style = f' style="background-color: {total_bg_color}; color: white; font-weight: bold;"' if total_bg_color else ' style="font-weight: bold;"'
            genome_list_html += f'<td class="count-cell"{total_style}>{total}</td>'

            for bgc_type in bgc_type_columns:
                count = bgc_types.get(bgc_type, 0)
                font_color = get_red_font_color(count, global_bgc_max)
                style = f' style="color: {font_color}; font-weight: bold;"' if font_color else ''
                genome_list_html += f'<td class="count-cell"{style}>{count if count > 0 else ""}</td>'

            genome_list_html += '</tr>'

        genome_list_html += '</table></div>'
        return genome_list_html

    def render_node(node, level=0):
        """Recursively render tree nodes."""
        name = node.get('name', 'Unknown')
        rank = node.get('rank', '')
        stats = node.get('stats', {})
        children = node.get('children', {})
        genomes = node.get('genomes', [])

        genome_count = stats.get('genome_count', 0)
        total_bgcs = stats.get('total_bgcs', 0)
        avg_bgcs = stats.get('avg_bgcs_per_genome', 0)
        genomes_with_bgcs = stats.get('genomes_with_bgcs', 0)
        bgc_distribution = stats.get('bgc_type_distribution', {})

        bgc_types_html = ''
        if bgc_distribution:
            bgc_types_list = []
            for bgc_type, type_stats in list(bgc_distribution.items())[:3]:
                count = type_stats.get('total_count', 0)
                pct = type_stats.get('percentage_of_genomes', 0)
                bgc_types_list.append(f"{bgc_type}: {count} ({pct:.0f}%)")
            bgc_types_html = ', '.join(bgc_types_list)
            if len(bgc_distribution) > 3:
                bgc_types_html += f" (+{len(bgc_distribution)-3} more)"

        node_id = f"node_{abs(hash(name + rank + str(level)))}"
        has_children = len(children) > 0
        has_genomes = len(genomes) > 0 and rank == 'species'

        toggle_icon = '&#9660;' if (has_children or has_genomes) else '&#9679;'
        toggle_class = 'toggle' if (has_children or has_genomes) else 'leaf'

        html = f'''
        <div class="tree-node level-{level}">
            <div class="node-header {toggle_class}" onclick="toggleNode('{node_id}')">
                <span class="toggle-icon">{toggle_icon}</span>
                <span class="node-name"><strong>{name}</strong> <em>({rank})</em></span>
                <span class="node-stats">
                    {genome_count} genomes | {total_bgcs} BGCs | avg {avg_bgcs:.2f}
                    {(' | ' + bgc_types_html) if bgc_types_html else ''}
                </span>
            </div>'''

        if has_children or has_genomes:
            html += f'<div class="node-children" id="{node_id}">'

            for child_name, child_node in children.items():
                html += render_node(child_node, level + 1)

            if has_genomes:
                html += render_genome_list(genomes, level + 1)

            html += '</div>'

        html += '</div>'
        return html

    tree_html = f'''
    <div class="taxonomy-tree-container">
        <style>
            .taxonomy-tree-container {{
                background: white;
                padding: 15px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 15px 0;
                font-size: 0.85em;
            }}
            .tree-node {{
                margin-left: 12px;  /* Consistent indent since nodes are nested */
            }}
            .tree-node.level-0 {{ margin-left: 0px; }}  /* Root has no indent */
            .node-header {{
                padding: 4px 6px;
                margin: 1px 0;
                cursor: pointer;
                border-radius: 3px;
                transition: background 0.2s;
                line-height: 1.3;
            }}
            .node-header:hover {{
                background: #f0f0f0;
            }}
            .node-header.leaf {{
                cursor: default;
            }}
            .toggle-icon {{
                display: inline-block;
                width: 14px;
                font-size: 0.8em;
                color: #2c5aa0;
            }}
            .node-name {{
                margin-right: 10px;
            }}
            .node-stats {{
                color: #666;
                font-size: 0.85em;
            }}
            .node-children {{
                margin-top: 1px;
            }}
            .genome-list {{
                margin: 6px 0 6px 16px;
                padding: 8px;
                background: #f8f9fa;
                border-radius: 4px;
                border-left: 2px solid #5b8ac5;
                overflow-x: auto;
                max-width: calc(100% - 20px);
            }}
            .genome-table {{
                width: 100%;
                border-collapse: collapse;
                font-size: 0.8em;
                background: white;
            }}
            .genome-table th {{
                background: #2c5aa0;
                color: white;
                padding: 4px 6px;
                text-align: left;
                font-weight: bold;
                border: 1px solid #ddd;
            }}
            .genome-table td {{
                padding: 4px 6px;
                border: 1px solid #ddd;
                text-align: center;
                transition: opacity 0.2s;
            }}
            .genome-table td:first-child {{
                text-align: left;
                background: white !important;
            }}
            .genome-table tr:hover td {{
                opacity: 0.8;
            }}
            .genome-link {{
                color: #2c5aa0;
                text-decoration: none;
                font-weight: 500;
            }}
            .genome-link:hover {{
                color: #5b8ac5;
                text-decoration: underline;
            }}
            .count-cell {{
                text-align: center;
                color: #999;
            }}
        </style>
        <script>
            function toggleNode(nodeId) {{
                const element = document.getElementById(nodeId);
                if (element) {{
                    const header = element.previousElementSibling;
                    const icon = header.querySelector('.toggle-icon');
                    if (element.style.display === 'none') {{
                        element.style.display = 'block';
                        icon.innerHTML = '&#9660;';
                    }} else {{
                        element.style.display = 'none';
                        icon.innerHTML = '&#9654;';
                    }}
                }}
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                const nodeChildren = document.querySelectorAll('.node-children');
                nodeChildren.forEach((node, index) => {{
                    if (index > 0) {{
                        node.style.display = 'none';
                        const icon = node.previousElementSibling.querySelector('.toggle-icon');
                        if (icon) icon.innerHTML = '&#9654;';
                    }}
                }});
            }});
        </script>
        {render_node(tree)}
    </div>'''

    return tree_html
