"""Visualization modules for BGC analysis pipeline."""

from .charts import get_bgc_color, plot_bgc_donut_chart, plot_kcb_identification_chart
from .tree_viz import (
    parse_newick,
    collect_phylo_leaves,
    prune_tree_to_leaves,
    tree_to_newick,
    calculate_phylo_positions,
    plot_circular_taxonomy_tree,
    plot_circular_phylogenetic_tree,
    prepare_phylo_tree_for_js,
    generate_static_circular_tree
)
from .tables import (
    get_genome_count,
    generate_genome_table_html,
    calculate_summary_statistics,
    create_bgc_distribution_table
)
from .clustering import (
    generate_bigslice_stats_html,
    generate_bigscape_stats_html,
    generate_gcf_visualization_html
)
from .taxonomy import generate_taxonomy_tree_html
from .resources import (
    aggregate_trace_by_process,
    generate_gantt_chart_html,
    generate_resource_usage_html
)

__all__ = [
    # charts
    'get_bgc_color',
    'plot_bgc_donut_chart',
    'plot_kcb_identification_chart',
    # tree_viz
    'parse_newick',
    'collect_phylo_leaves',
    'prune_tree_to_leaves',
    'tree_to_newick',
    'calculate_phylo_positions',
    'plot_circular_taxonomy_tree',
    'plot_circular_phylogenetic_tree',
    'prepare_phylo_tree_for_js',
    'generate_static_circular_tree',
    # tables
    'get_genome_count',
    'generate_genome_table_html',
    'calculate_summary_statistics',
    'create_bgc_distribution_table',
    # clustering
    'generate_bigslice_stats_html',
    'generate_bigscape_stats_html',
    'generate_gcf_visualization_html',
    # taxonomy
    'generate_taxonomy_tree_html',
    # resources
    'aggregate_trace_by_process',
    'generate_gantt_chart_html',
    'generate_resource_usage_html',
]
