/**
 * Download antiSMASH databases (~50 GB)
 * Required for BGC detection with antiSMASH
 */
process DOWNLOAD_ANTISMASH_DBS {
    tag "antiSMASH databases"
    label 'process_low'
    storeDir "${params.outdir}/databases"

    output:
    path "antismash/", emit: db_dir

    script:
    """
    echo "=============================================="
    echo "Downloading antiSMASH Databases"
    echo "=============================================="
    echo "Size: ~50 GB"
    echo "Started at: \$(date)"
    echo ""

    # Set the database directory
    export ANTISMASH_DB_PATH=\$(pwd)/antismash
    mkdir -p \$ANTISMASH_DB_PATH

    # Download databases
    echo "Downloading databases..."
    download-antismash-databases --database-dir \$ANTISMASH_DB_PATH

    # Verify download
    echo ""
    echo "Verifying download..."
    if [ ! -d "\$ANTISMASH_DB_PATH" ] || [ -z "\$(ls -A \$ANTISMASH_DB_PATH)" ]; then
        echo "ERROR: Database directory is empty"
        exit 1
    fi

    echo ""
    echo "Database contents:"
    ls -lh \$ANTISMASH_DB_PATH/

    echo ""
    echo "=============================================="
    echo "antiSMASH databases ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
