process TABULATE_REGIONS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path "antismash_results/*"

    output:
    path "region_tabulation.tsv", emit: tabulation

    script:
    """
    python ${projectDir}/scripts/analysis/tabulate_regions.py antismash_results region_tabulation.tsv
    """
}
