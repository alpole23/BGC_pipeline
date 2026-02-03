/**
 * Check if antiSMASH results can be reused from a previous taxon run.
 * Validates that:
 * 1. Results exist for the genome
 * 2. The antiSMASH version matches (major.minor)
 * 3. The parameter configuration matches
 */
process CHECK_ANTISMASH_REUSE {
    tag "$genome.baseName"
    label 'process_local'
    executor 'local'
    cache false  // Always check filesystem

    input:
    val taxon
    val reuse_taxon
    path genome
    val antismash_version
    val antismash_params_hash

    output:
    tuple path(genome), env(STATUS), env(EXISTING_PATH), emit: check_result

    script:
    def genome_name = genome.baseName
    // Use Utils helper to build absolute reuse path
    def reuse_dir = Utils.buildReusePath(params, projectDir, "antismash", reuse_taxon, genome_name)
    def meta_file = "${reuse_dir}/.antismash_meta"
    """
    STATUS="RUN"
    EXISTING_PATH=""

    if [ -d "${reuse_dir}" ] && [ -f "${reuse_dir}/${genome_name}.json" ]; then
        if [ -f "${meta_file}" ]; then
            # Read stored version and params hash
            stored_version=\$(grep "^version=" "${meta_file}" | cut -d= -f2)
            stored_params=\$(grep "^params_hash=" "${meta_file}" | cut -d= -f2)

            # Compare major.minor version (ignore patch)
            stored_mm=\$(echo "\$stored_version" | cut -d. -f1,2)
            current_mm=\$(echo "${antismash_version}" | cut -d. -f1,2)

            if [ "\$stored_mm" = "\$current_mm" ] && [ "\$stored_params" = "${antismash_params_hash}" ]; then
                STATUS="REUSE"
                EXISTING_PATH="${reuse_dir}"
                echo "REUSE: ${genome_name} - version \$stored_mm matches, params match"
            else
                echo "RE-RUN: ${genome_name} - version or params mismatch (stored=\$stored_version/\$stored_params current=${antismash_version}/${antismash_params_hash})"
            fi
        else
            # No meta file = legacy result, still reuse but warn
            STATUS="REUSE"
            EXISTING_PATH="${reuse_dir}"
            echo "REUSE (legacy): ${genome_name} - no meta file, assuming compatible"
        fi
    else
        echo "RUN: ${genome_name} - no existing results found"
    fi
    """
}

/**
 * Copy antiSMASH results from a previous taxon directory to the current one.
 */
process COPY_ANTISMASH_RESULT {
    tag "$genome_name"
    label 'process_low'
    publishDir "${params.outdir}/antismash_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon
    val genome_name
    val existing_result_path  // Use val instead of path to avoid symlink staging issues

    output:
    path "${genome_name}/", emit: result_dir

    script:
    """
    cp -rL "${existing_result_path}" "${genome_name}"
    """
}
