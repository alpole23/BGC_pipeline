process COUNT_REGIONS {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path "antismash_results/*"

    output:
    path "region_counts.tsv", emit: counts

    script:
    def by_contig = params.count_per_contig ? "--by_contig" : ""
    def split_hybrids = params.split_hybrids ? "--split_hybrids" : ""
    """
    python ${projectDir}/scripts/analysis/count_regions.py antismash_results region_counts.tsv ${by_contig} ${split_hybrids}
    """
}
