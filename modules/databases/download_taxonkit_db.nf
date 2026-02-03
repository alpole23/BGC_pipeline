/**
 * Download NCBI taxonomy database (~100 MB)
 * Required for TaxonKit taxonomy processing
 */
process DOWNLOAD_TAXONKIT_DB {
    tag "NCBI taxonomy"
    label 'process_low'
    storeDir "${params.outdir}/databases/taxonkit"

    output:
    path "taxdump", emit: taxdump_dir

    script:
    """
    echo "=============================================="
    echo "Downloading NCBI Taxonomy Database"
    echo "=============================================="
    echo "Size: ~100 MB"
    echo "Started at: \$(date)"
    echo ""

    # Create directory for taxonomy dump
    mkdir -p taxdump

    # Download NCBI taxonomy dump
    echo "Downloading taxdump.tar.gz from NCBI FTP..."
    wget -q --show-progress -O taxdump.tar.gz "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

    # Extract the archive
    echo ""
    echo "Extracting taxonomy files..."
    tar -xzf taxdump.tar.gz -C taxdump

    # Clean up archive
    rm taxdump.tar.gz

    # Verify required files
    echo ""
    echo "Verifying download..."
    REQUIRED_FILES="nodes.dmp names.dmp delnodes.dmp merged.dmp"
    for f in \$REQUIRED_FILES; do
        if [ ! -f "taxdump/\$f" ]; then
            echo "ERROR: Required file taxdump/\$f not found"
            exit 1
        fi
    done

    echo ""
    echo "Database contents:"
    ls -lh taxdump/

    echo ""
    echo "=============================================="
    echo "NCBI taxonomy database ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
