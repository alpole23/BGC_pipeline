process VISUALIZE_RESULTS {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path counts_file
    path tabulation_file
    path assembly_info
    path name_map
    path taxonomy_map
    path taxonomy_tree
    path bigscape_stats
    path bigscape_db
    path gcf_data
    path phylo_tree
    path gtdbtk_summary
    path trace_file
    path versions_file
    path gcf_tree_png

    output:
    path "*.png", emit: plots, optional: true
    path "*.html", emit: reports, optional: true
    path "*.nwk", emit: newick_files, optional: true
    path "genomes/*.html", emit: genome_pages, optional: true

    script:

    def counts_arg = counts_file.name != 'NO_FILE' ? "--counts ${counts_file}" : ""
    def tab_arg = tabulation_file.name != 'NO_FILE' ? "--tabulation ${tabulation_file}" : ""
    def assembly_arg = assembly_info.name != 'NO_FILE' ? "--assembly_info ${assembly_info}" : ""
    def name_map_arg = name_map.name != 'NO_FILE' ? "--name_map ${name_map}" : ""
    def taxonomy_map_arg = taxonomy_map && taxonomy_map.name != 'NO_FILE' ? "--taxonomy_map ${taxonomy_map}" : ""
    def taxonomy_tree_arg = taxonomy_tree && taxonomy_tree.name != 'NO_FILE' ? "--taxonomy_tree ${taxonomy_tree}" : ""
    def bigscape_stats_arg = Utils.isValidInput(bigscape_stats) ? "--bigscape_stats ${Utils.getFirstFile(bigscape_stats)}" : ""
    def bigscape_db_arg = Utils.isValidInput(bigscape_db) ? "--bigscape_db ${Utils.getFirstFile(bigscape_db)}" : ""
    def gcf_data_arg = Utils.isValidInput(gcf_data) ? "--gcf_data ${Utils.getFirstFile(gcf_data)}" : ""
    def phylo_tree_arg = Utils.isValidInput(phylo_tree) ? "--phylo_tree ${Utils.getFirstFile(phylo_tree)}" : ""
    def gtdbtk_summary_arg = Utils.isValidInput(gtdbtk_summary) ? "--gtdbtk_summary ${Utils.getFirstFile(gtdbtk_summary)}" : ""
    def trace_arg = Utils.isValidInput(trace_file) ? "--trace ${Utils.getFirstFile(trace_file)}" : ""
    def versions_arg = Utils.isValidInput(versions_file) ? "--versions ${Utils.getFirstFile(versions_file)}" : ""
    def mibig_arg = params.bigscape_mibig_version ? "--mibig_included" : ""
    def skip_tree_arg = params.skip_tree ? "--skip_tree" : ""
    def outgroup_arg = params.gtdbtk_outgroup ? "--outgroup '${params.gtdbtk_outgroup}'" : ""
    def gcf_tree_arg = gcf_tree_png.name != 'NO_GCF_TREE' ? "--gcf_tree ${gcf_tree_png}" : ""
    """
    python ${projectDir}/scripts/visualize_results.py ${counts_arg} ${tab_arg} ${assembly_arg} ${name_map_arg} ${taxonomy_map_arg} ${taxonomy_tree_arg} ${bigscape_stats_arg} ${bigscape_db_arg} ${gcf_data_arg} ${phylo_tree_arg} ${gtdbtk_summary_arg} ${trace_arg} ${versions_arg} ${mibig_arg} ${skip_tree_arg} ${outgroup_arg} ${gcf_tree_arg} --outdir . --taxon "${taxon}"
    """
}
