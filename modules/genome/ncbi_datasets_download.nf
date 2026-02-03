process NCBI_DATASETS_DOWNLOAD {
    tag "$taxon"
    label 'process_medium'
    label 'retry_on_error'
    publishDir "${params.outdir}/ncbi_genomes/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'

    input:
    val taxon

    output:
    path "ncbi_dataset/data/**/*.gbff", emit: genomes
    path "ncbi_dataset/", emit: dataset_dir
    path "ncbi_dataset/data/assembly_info_table.txt", emit: assembly_info
    path "ncbi_dataset/data/assembly_data_report.jsonl", emit: assembly_data_report
    path "ncbi_dataset/data/taxonomy_report.jsonl", emit: taxonomy_report, optional: true

    script:
    """    
    datasets download genome taxon "${taxon}" \
        --include gbff \
        --assembly-source GenBank \
        --exclude-atypical \
        --filename ncbi_dataset.zip \
        --dehydrated

    unzip -q ncbi_dataset.zip

    # Rehydrate with retry logic for transient network errors
    MAX_RETRIES=3
    RETRY_COUNT=0
    until datasets rehydrate --directory . ; do
        RETRY_COUNT=\$((RETRY_COUNT + 1))
        if [ \$RETRY_COUNT -ge \$MAX_RETRIES ]; then
            echo "ERROR: Rehydrate failed after \$MAX_RETRIES attempts"
            exit 1
        fi
        echo "Rehydrate attempt \$RETRY_COUNT failed, retrying in 30 seconds..."
        sleep 30
    done

    # Validate downloaded files and re-rehydrate corrupted ones
    # Corrupted files from failed downloads contain only null bytes
    # Force filesystem sync to ensure all writes are flushed before validation
    echo "Syncing filesystem before validation..."
    sync
    sleep 5

    echo "Validating downloaded genome files..."
    MAX_VALIDATION_CYCLES=3
    VALIDATION_CYCLE=0

    while [ \$VALIDATION_CYCLE -lt \$MAX_VALIDATION_CYCLES ]; do
        VALIDATION_CYCLE=\$((VALIDATION_CYCLE + 1))
        CORRUPTED_COUNT=0
        TOTAL_COUNT=0

        # Check each .gbff file for corruption (first 100 bytes all null = corrupted)
        # Use cat to force actual disk read, bypassing any page cache issues
        for gbff in ncbi_dataset/data/*/*.gbff; do
            if [ -f "\$gbff" ]; then
                TOTAL_COUNT=\$((TOTAL_COUNT + 1))
                # Drop file from page cache and re-read to get actual disk content
                # Valid GenBank files start with "LOCUS" - check for any printable chars
                FIRST_BYTES=\$(cat "\$gbff" 2>/dev/null | head -c 100 | tr -d '\0' | head -c 5)
                if [ -z "\$FIRST_BYTES" ]; then
                    echo "Corrupted file detected: \$gbff"
                    # Delete the corrupted file so rehydrate will re-download it
                    rm -f "\$gbff"
                    CORRUPTED_COUNT=\$((CORRUPTED_COUNT + 1))
                fi
            fi
        done

        echo "Validation cycle \$VALIDATION_CYCLE: \$CORRUPTED_COUNT/\$TOTAL_COUNT files corrupted"

        if [ \$CORRUPTED_COUNT -eq 0 ]; then
            echo "All files validated successfully"
            break
        fi

        if [ \$VALIDATION_CYCLE -lt \$MAX_VALIDATION_CYCLES ]; then
            echo "Re-rehydrating \$CORRUPTED_COUNT corrupted files..."
            sleep 10
            datasets rehydrate --directory . || {
                echo "WARNING: Re-rehydration attempt \$VALIDATION_CYCLE failed"
            }
            # Sync after re-rehydration to ensure writes are flushed
            sync
            sleep 5
        else
            echo "WARNING: \$CORRUPTED_COUNT files still corrupted after \$MAX_VALIDATION_CYCLES validation cycles"
            echo "These genomes will be skipped in downstream analysis"
        fi
    done

    # Final count of valid files (use cat to force disk read)
    VALID_COUNT=0
    for gbff in ncbi_dataset/data/*/*.gbff; do
        if [ -f "\$gbff" ]; then
            FIRST_BYTES=\$(cat "\$gbff" 2>/dev/null | head -c 5 | tr -d '\0')
            if [ -n "\$FIRST_BYTES" ]; then
                VALID_COUNT=\$((VALID_COUNT + 1))
            fi
        fi
    done
    TOTAL_FILES=\$(find ncbi_dataset/data -name "*.gbff" -type f | wc -l)
    echo "Final count: \$VALID_COUNT/\$TOTAL_FILES valid genome files"

    if [ \$VALID_COUNT -eq 0 ]; then
        echo "ERROR: No valid genome files after rehydration"
        exit 1
    fi

    dataformat tsv genome \
        --inputfile ncbi_dataset/data/assembly_data_report.jsonl \
        --fields accession,organism-name,organism-infraspecific-strain,assminfo-biosample-isolation-source,assminfo-biosample-isolate,assminfo-notes,checkm-completeness,checkm-contamination,organism-tax-id \
        > ncbi_dataset/data/assembly_info_table.txt

    # Download taxonomy data using the user-specified taxon
    echo "Downloading taxonomy data for taxon: ${taxon}"
    datasets download taxonomy taxon "${taxon}" \\
        --filename taxonomy.zip \\
        2>/dev/null || {
        echo "WARNING: Taxonomy download failed for taxon ${taxon}"
        exit 0
    }

    # Unzip taxonomy data
    unzip -o -q taxonomy.zip

    # Copy the pre-generated TSV file if available
    if [ -f "ncbi_dataset/data/taxonomy_summary.tsv" ]; then
        cp ncbi_dataset/data/taxonomy_summary.tsv taxonomy_report.tsv
        echo "Taxonomy report generated successfully"
    else
        echo "WARNING: taxonomy_summary.tsv not found in downloaded archive"
    fi
    """

}