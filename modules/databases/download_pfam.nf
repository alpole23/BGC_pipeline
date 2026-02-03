/**
 * Download Pfam-A HMM database (~1.5 GB)
 * Required for BiG-SCAPE clustering
 */
process DOWNLOAD_PFAM {
    tag "Pfam database"
    label 'process_low'
    storeDir "${params.outdir}/databases"

    output:
    path "pfam/", emit: pfam_db

    script:
    """
    echo "=============================================="
    echo "Downloading Pfam Database"
    echo "=============================================="
    echo "Size: ~1.5 GB"
    echo "Started at: \$(date)"
    echo ""

    mkdir -p pfam
    cd pfam

    # Download Pfam-A HMM profiles
    if [ ! -f "Pfam-A.hmm.gz" ] && [ ! -f "Pfam-A.hmm" ]; then
        echo "Downloading Pfam-A.hmm..."
        wget -q --show-progress https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        echo "Extracting..."
        gunzip Pfam-A.hmm.gz
    else
        echo "Pfam-A.hmm already exists"
    fi

    # Press the HMM database for faster searching
    if [ ! -f "Pfam-A.hmm.h3p" ]; then
        echo "Pressing Pfam database with hmmpress..."
        hmmpress Pfam-A.hmm
    else
        echo "Pfam database already pressed"
    fi

    # Verify
    echo ""
    echo "Verifying download..."
    if [ ! -f "Pfam-A.hmm" ] || [ ! -f "Pfam-A.hmm.h3p" ]; then
        echo "ERROR: Pfam database files missing"
        exit 1
    fi

    echo ""
    echo "Database contents:"
    ls -lh

    echo ""
    echo "=============================================="
    echo "Pfam database ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
