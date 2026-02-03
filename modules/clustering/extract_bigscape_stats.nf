process EXTRACT_BIGSCAPE_STATS {
    tag "$taxon"
    conda "conda-forge::python=3.9"
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path bigscape_dir

    output:
    path "bigscape_statistics.json", emit: stats_json

    script:
    """
    python ${projectDir}/scripts/clustering/extract_bigscape_stats.py ${bigscape_dir} bigscape_statistics.json
    """
}
