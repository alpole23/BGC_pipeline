process EXTRACT_BIGSLICE_STATS {
    tag "$taxon"
    conda "conda-forge::sqlite conda-forge::python=3.9"
    publishDir "${params.outdir}/bigslice_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path bigslice_dir

    output:
    path "bigslice_statistics.json", emit: stats_json

    script:
    """
    python ${projectDir}/scripts/clustering/extract_bigslice_stats.py ${bigslice_dir} bigslice_statistics.json
    """
}
