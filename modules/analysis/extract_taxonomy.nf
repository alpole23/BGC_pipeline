process EXTRACT_TAXONOMY {
    tag "Extracting taxonomy from NCBI metadata"
    label 'process_low'
    publishDir "${params.outdir}/main_analysis_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path assembly_report
    path taxonomy_report
    path taxdump_dir

    output:
    path "taxonomy_map.json", emit: taxonomy_map

    script:
    """
    export TAXONKIT_DB='${taxdump_dir}'
    python ${projectDir}/scripts/taxonomy/extract_taxonomy.py ${assembly_report} ${taxdump_dir} taxonomy_map.json
    """
}
