import java.security.MessageDigest

/**
 * Shared utility functions for the BGC pipeline.
 * These are automatically loaded by Nextflow from the lib/ directory.
 *
 * Usage in modules: Utils.sanitizeTaxon(name)
 * Usage in main.nf: sanitizeTaxon(name) - uses local function definition
 */
class Utils {
    /**
     * Sanitize taxon name for use in file paths and directory names.
     * Replaces special characters with underscores for filesystem safety.
     */
    static String sanitizeTaxon(String name) {
        name.replaceAll('[^a-zA-Z0-9_]', '_').replaceAll('_+', '_').replaceAll('^_|_$', '')
    }

    /**
     * Generate MD5 hash of relevant antiSMASH parameters.
     * Used to detect when parameters change and results need recomputation.
     * Note: clusterhmmer and tigrfam are always enabled, not tracked.
     *
     * @param params The Nextflow params object
     * @return MD5 hash string of the parameter configuration
     */
    static String antismashParamsHash(params) {
        // Collect relevant antiSMASH parameters into a sorted string
        def relevantParams = [
            "minimal=${params.antismash_minimal ?: false}",
            "cb_general=${params.antismash_cb_general ?: false}",
            "cc_mibig=${params.antismash_cc_mibig ?: false}",
            "cb_knownclusters=${params.antismash_cb_knownclusters ?: true}",
            "smcog_trees=${params.antismash_smcog_trees ?: false}",
            "hmmdetection_rules=${params.hmmdetection_rules ?: ''}"
        ].sort().join(";")

        // Return MD5 hash of the parameter string
        return md5(relevantParams)
    }

    /**
     * Calculate MD5 hash of a string.
     */
    static String md5(String input) {
        MessageDigest md = MessageDigest.getInstance("MD5")
        byte[] digest = md.digest(input.bytes)
        return digest.collect { String.format("%02x", it) }.join()
    }

    /**
     * Build absolute path for reusing results from a previous taxon run.
     * Handles both absolute and relative outdir paths.
     *
     * @param params The Nextflow params object (needs outdir)
     * @param projectDir The Nextflow projectDir variable
     * @param toolName The tool name (e.g., "antismash", "gtdbtk")
     * @param reuseTaxon The taxon name to reuse results from
     * @param subPath Optional sub-path within the result directory
     * @return Absolute path to the reuse directory
     */
    static String buildReusePath(params, projectDir, String toolName, String reuseTaxon, String subPath = "") {
        def outdir = params.outdir?.toString() ?: "results"
        def outdir_abs = outdir.startsWith('/') ? outdir : "${projectDir}/${outdir}"
        def basePath = "${outdir_abs}/${toolName}_results/${sanitizeTaxon(reuseTaxon)}"
        return subPath ? "${basePath}/${subPath}" : basePath
    }

    /**
     * Check if a file path represents a valid input (not a placeholder).
     * Handles both single files and lists from glob patterns.
     *
     * @param input The file path or list of paths to check
     * @return true if the input is valid, false if it's a placeholder
     */
    static boolean isValidInput(input) {
        if (!input) return false
        if (input instanceof List) {
            return input.size() > 0 && !input[0].name?.startsWith('NO_')
        }
        return !input.name?.startsWith('NO_')
    }
}
