process BIGSCAPE {
    tag "$taxon"
    label 'process_high'
    publishDir "${params.outdir}/bigscape_results", mode: 'copy'
    cache 'lenient'  // Use lenient caching for directory inputs

    input:
    val taxon
    path antismash_results
    path pfam_db
    
    output:
    path "${Utils.sanitizeTaxon(taxon)}/", emit: bigscape_dir
    path "${Utils.sanitizeTaxon(taxon)}/${Utils.sanitizeTaxon(taxon)}.db", emit: bigscape_db

    script:
    def safe_taxon = Utils.sanitizeTaxon(taxon)
    def cutoffs = params.bigscape_cutoffs ?: "0.30"
    def alignment_mode = params.bigscape_alignment_mode ?: "auto"
    def mibig_version = params.bigscape_mibig_version ? "--mibig-version ${params.bigscape_mibig_version}" : ""
    def classify = params.bigscape_classify ? "--classify ${params.bigscape_classify}" : ""
    def include_singletons = params.bigscape_include_singletons ? "--include-singletons" : ""
    def mix = params.bigscape_mix ? "--mix" : ""
    """
    # BiG-SCAPE 2 requires path to Pfam-A.hmm file
    export PFAM_PATH=\$(readlink -f ${pfam_db})/Pfam-A.hmm
    
    echo "Running BiG-SCAPE on ${taxon}"
    echo "Using Pfam database: \$PFAM_PATH"
    echo "Using ${task.cpus} cores"
    
    # Create input directory with all antiSMASH results
    mkdir -p antismash_input
    
    # Copy/link all antiSMASH result directories
    for dir in ${antismash_results}; do
        if [ -d "\$dir" ]; then
            echo "Linking \$dir"
            ln -s \$(readlink -f "\$dir") antismash_input/
        fi
    done
    
    # Count how many we have
    INPUT_COUNT=\$(ls -d antismash_input/*/ 2>/dev/null | wc -l)
    echo "Processing \$INPUT_COUNT antiSMASH results"
    
    if [ "\$INPUT_COUNT" -eq 0 ]; then
        echo "ERROR: No antiSMASH results found"
        echo "Listing what we received:"
        ls -la
        exit 1
    fi
    
    # Run BiG-SCAPE 2
    bigscape cluster \
        -i antismash_input \
        -o ${safe_taxon} \
        --pfam-path \$PFAM_PATH \
        --alignment-mode ${alignment_mode} \
        --gcf-cutoffs ${cutoffs} \
        ${include_singletons} \
        ${mix} \
        ${mibig_version} \
        --cores ${task.cpus} \
        ${classify}

    # Update index.html to use taxon-specific database name
    # BiG-SCAPE generates index.html with hardcoded 'data_sqlite.db' reference
    # Replace it with the actual taxon-specific database filename
    cd ${safe_taxon}
    sed -i "s|data_sqlite.db|${safe_taxon}.db|g" index.html
    cd ..

    echo "BiG-SCAPE analysis complete"
    """
}