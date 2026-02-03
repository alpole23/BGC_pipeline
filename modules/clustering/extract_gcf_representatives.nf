process EXTRACT_GCF_REPRESENTATIVES {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/bigscape_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path bigscape_dir
    path "antismash_input/*", stageAs: 'antismash_input/*'
    path tabulation_file

    output:
    path "gcf_representatives.json", emit: gcf_data

    script:
    def taxon_clean = Utils.sanitizeTaxon(taxon)
    def tabulation_arg = tabulation_file.name != 'NO_TABULATION' ? "--tabulation ${tabulation_file}" : ""
    """
    python ${projectDir}/scripts/clustering/extract_gcf_representatives.py ${bigscape_dir} antismash_input gcf_representatives.json --taxon "${taxon_clean}" ${tabulation_arg}
    """
}
