process GCF_BIOSYNTHETIC_TREE {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}/gcf_heatmap", mode: 'copy'

    input:
    val taxon
    path bigscape_db
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path gtdbtk_tree
    path gtdbtk_summary

    output:
    path "gcf_species_heatmap.png",       emit: heatmap_png,   optional: true
    path "gcf_species_heatmap.svg",       emit: heatmap_svg,   optional: true
    path "gcf_biosynthetic_tree.png",     emit: gcf_tree_png,  optional: true
    path "gcf_biosynthetic_tree.svg",     emit: gcf_tree_svg,  optional: true
    path "phosphonate_metadata.json",     emit: metadata,      optional: true

    script:
    def tree_arg    = gtdbtk_tree.name    != 'NO_PHYLO_TREE'     ? "--gtdbtk_tree ${gtdbtk_tree}"      : ""
    def summary_arg = gtdbtk_summary.name != 'NO_GTDBTK_SUMMARY' ? "--gtdbtk_summary ${gtdbtk_summary}" : ""
    """
    # Step 1: Generate BGC metadata from BiG-SCAPE database (metadata only — skip slow NJ tree)
    python ${projectDir}/scripts/bgc_pfam_tree.py \\
        --db ${bigscape_db} \\
        --bgc_type phosphonate \\
        --metadata_only \\
        --outdir .

    # Step 2: Generate coupling enzyme class annotations
    python ${projectDir}/scripts/bgc_coupling_annotation.py \\
        --antismash_dir antismash_input \\
        --metadata phosphonate_metadata.json \\
        --outfile phosphonate_itol_coupling.txt \\
        --bgc_type phosphonate

    # Step 3: Generate GCF × species heatmap
    python ${projectDir}/scripts/bgc_gcf_heatmap.py \\
        --metadata phosphonate_metadata.json \\
        --coupling_annotation phosphonate_itol_coupling.txt \\
        ${tree_arg} \\
        ${summary_arg} \\
        --outdir .

    # Step 4: Generate GCF biosynthetic NJ tree figure
    python ${projectDir}/scripts/bgc_gcf_tree.py \\
        --db ${bigscape_db} \\
        --coupling_annotation phosphonate_itol_coupling.txt \\
        --outdir .
    """
}
