process RENAME_GENOMES {
    tag "$assembly_id"
    label 'process_low'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(taxon)}/renamed_genomes", mode: 'copy'

    input:
    val taxon
    tuple val(assembly_id), path(genome_file)
    path name_map

    output:
    path "*.gbff", emit: renamed_genome

    script:
    """
    python ${projectDir}/scripts/genome/rename_genome.py "${assembly_id}" ${genome_file} ${name_map}
    """
}
