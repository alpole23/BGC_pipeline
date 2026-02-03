process CREATE_NAME_MAP {
    tag "$taxon"
    label 'process_low'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    path assembly_info

    output:
    path "name_map.json", emit: name_map

    script:
    """
    python ${projectDir}/scripts/genome/create_name_map.py ${assembly_info} name_map.json
    """
}
