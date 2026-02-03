process VISUALIZE_RESULTS {
    tag "$taxon"
    label 'process_medium'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(taxon)}/main_data_visualization", mode: 'copy'

    input:
    val taxon
    path counts_file
    path tabulation_file
    path assembly_info
    path name_map
    path taxonomy_map
    path taxonomy_tree
    path bigslice_stats
    path bigscape_stats
    path bigscape_db
    path gcf_data
    path phylo_tree
    path gtdbtk_summary
    path trace_file
    path versions_file

    output:
    path "*.png", emit: plots, optional: true
    path "*.html", emit: reports, optional: true
    path "*.nwk", emit: newick_files, optional: true
    path "genomes/*.html", emit: genome_pages, optional: true

    script:
    // Helper to check if input is valid (not a placeholder)
    // Handles both single files and lists (from glob patterns)
    def isValidInput = { input ->
        if (!input) return false
        if (input instanceof List) {
            return input.size() > 0 && !input[0].name.startsWith('NO_')
        }
        return !input.name.startsWith('NO_')
    }
    // Helper to get first file from input (handles lists)
    def getFirstFile = { input ->
        if (input instanceof List) return input[0]
        return input
    }

    def counts_arg = counts_file.name != 'NO_FILE' ? "--counts ${counts_file}" : ""
    def tab_arg = tabulation_file.name != 'NO_FILE' ? "--tabulation ${tabulation_file}" : ""
    def assembly_arg = assembly_info.name != 'NO_FILE' ? "--assembly_info ${assembly_info}" : ""
    def name_map_arg = name_map.name != 'NO_FILE' ? "--name_map ${name_map}" : ""
    def taxonomy_map_arg = taxonomy_map && taxonomy_map.name != 'NO_FILE' ? "--taxonomy_map ${taxonomy_map}" : ""
    def taxonomy_tree_arg = taxonomy_tree && taxonomy_tree.name != 'NO_FILE' ? "--taxonomy_tree ${taxonomy_tree}" : ""
    def bigslice_stats_arg = isValidInput(bigslice_stats) ? "--bigslice_stats ${getFirstFile(bigslice_stats)}" : ""
    def bigscape_stats_arg = isValidInput(bigscape_stats) ? "--bigscape_stats ${getFirstFile(bigscape_stats)}" : ""
    def bigscape_db_arg = isValidInput(bigscape_db) ? "--bigscape_db ${getFirstFile(bigscape_db)}" : ""
    def gcf_data_arg = isValidInput(gcf_data) ? "--gcf_data ${getFirstFile(gcf_data)}" : ""
    def phylo_tree_arg = isValidInput(phylo_tree) ? "--phylo_tree ${getFirstFile(phylo_tree)}" : ""
    def gtdbtk_summary_arg = isValidInput(gtdbtk_summary) ? "--gtdbtk_summary ${getFirstFile(gtdbtk_summary)}" : ""
    def trace_arg = isValidInput(trace_file) ? "--trace ${getFirstFile(trace_file)}" : ""
    def versions_arg = isValidInput(versions_file) ? "--versions ${getFirstFile(versions_file)}" : ""
    def mibig_arg = params.bigscape_mibig_version ? "--mibig_included" : ""
    def outgroup_arg = params.gtdbtk_outgroup ? "--outgroup '${params.gtdbtk_outgroup}'" : ""
    """
    python ${projectDir}/scripts/visualize_results.py ${counts_arg} ${tab_arg} ${assembly_arg} ${name_map_arg} ${taxonomy_map_arg} ${taxonomy_tree_arg} ${bigslice_stats_arg} ${bigscape_stats_arg} ${bigscape_db_arg} ${gcf_data_arg} ${phylo_tree_arg} ${gtdbtk_summary_arg} ${trace_arg} ${versions_arg} ${mibig_arg} ${outgroup_arg} --outdir . --taxon "${taxon}"
    """
}
