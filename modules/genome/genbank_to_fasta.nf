process GENBANK_TO_FASTA {
    tag "$genome.baseName"
    label 'process_low'
    label 'tolerant'

    input:
    path genome

    output:
    path "*.fna", emit: fasta

    script:
    def output_name = genome.baseName + ".fna"
    """
    python ${projectDir}/scripts/genome/genbank_to_fasta.py ${genome} ${output_name}
    """
}
