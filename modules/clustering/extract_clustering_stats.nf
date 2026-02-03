/**
 * Extract statistics from clustering tool output (BiG-SCAPE or BiG-SLiCE).
 * This is a parameterized process that handles both tools.
 */
process EXTRACT_CLUSTERING_STATS {
    tag "$taxon - $tool"
    label 'process_low'
    publishDir "${params.outdir}/${tool}_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    val tool       // "bigscape" or "bigslice"
    path input_dir

    output:
    path "${tool}_statistics.json", emit: stats_json

    script:
    """
    python ${projectDir}/scripts/clustering/extract_${tool}_stats.py ${input_dir} ${tool}_statistics.json
    """
}
