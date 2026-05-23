/**
 * Extract statistics from BiG-SCAPE clustering output.
 */
process EXTRACT_CLUSTERING_STATS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(params.taxon)}", mode: 'copy'

    input:
    val taxon
    path input_dir

    output:
    path "bigscape_statistics.json", emit: stats_json

    script:
    """
    python ${projectDir}/scripts/clustering/extract_bigscape_stats.py ${input_dir} bigscape_statistics.json
    """
}
