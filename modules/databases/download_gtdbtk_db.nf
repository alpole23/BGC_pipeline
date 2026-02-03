/**
 * Download GTDB-Tk reference database (~140 GB)
 * Required for phylogenetic classification with GTDB-Tk
 * Uses parallel downloads via aria2 for faster completion
 */
process DOWNLOAD_GTDBTK_DB {
    tag "GTDB-Tk database"
    label 'process_medium'
    storeDir "${params.outdir}/databases/gtdbtk"

    output:
    path "release*", emit: db_dir

    script:
    """
    echo "=============================================="
    echo "GTDB-Tk Database Download (Split Package - Parallel)"
    echo "=============================================="
    echo "WARNING: This database is approximately 140 GB"
    echo "Downloading 14 x 10GB split files in parallel for maximum speed"
    echo ""

    # AAU mirror split package - downloads all parts in parallel
    BASE_URL="https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package"

    echo "Downloading from: \$BASE_URL"
    echo "Started at: \$(date)"
    echo ""

    # Create download list for aria2 (all 14 parts)
    cat > download_list.txt << 'URLS'
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_aa
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ab
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ac
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ad
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ae
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_af
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ag
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ah
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ai
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_aj
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_ak
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_al
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_am
https://data.gtdb.aau.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/split_package/gtdbtk_data.tar.gz.part_an
URLS

    # Download all 14 parts in parallel
    # -j 14: download 14 files simultaneously
    # -x 4: 4 connections per file (14 files x 4 = 56 total connections)
    # -s 4: split each file into 4 segments
    aria2c -j 14 -x 4 -s 4 -k 20M \\
        --file-allocation=none \\
        --console-log-level=notice \\
        --summary-interval=30 \\
        -i download_list.txt

    echo ""
    echo "All parts downloaded at: \$(date)"
    echo "Part files:"
    ls -lh gtdbtk_data.tar.gz.part_*

    echo ""
    echo "Concatenating split files..."
    cat gtdbtk_data.tar.gz.part_* > gtdbtk_data.tar.gz

    echo "Combined file size: \$(du -h gtdbtk_data.tar.gz)"

    # Clean up split files
    rm -f gtdbtk_data.tar.gz.part_*
    rm -f download_list.txt

    echo ""
    echo "Extracting database (this may take 30-60 minutes)..."
    tar -xzf gtdbtk_data.tar.gz

    # Clean up archive to save space
    rm gtdbtk_data.tar.gz

    # Find the release directory (e.g., release220, release226)
    RELEASE_DIR=\$(ls -d release* 2>/dev/null | head -1)

    if [ -z "\$RELEASE_DIR" ]; then
        echo "ERROR: Could not find extracted GTDB-Tk release directory"
        exit 1
    fi

    echo ""
    echo "GTDB-Tk database downloaded successfully to \$RELEASE_DIR"
    echo "Database contents:"
    ls -lh \$RELEASE_DIR/

    echo ""
    echo "=============================================="
    echo "GTDB-Tk database ready for use"
    echo "=============================================="
    """
}
