/**
 * Download BiG-SLiCE HMM models (~500 MB)
 * Required for BiG-SLiCE clustering
 */
process DOWNLOAD_BIGSLICE_MODELS {
    tag "BiG-SLiCE models"
    label 'process_low'
    storeDir "${params.outdir}/databases/bigslice_models"

    output:
    path "hmmdb/", emit: models_dir

    script:
    """
    echo "=============================================="
    echo "Setting up BiG-SLiCE HMM Models"
    echo "=============================================="
    echo "Size: ~500 MB"
    echo "Started at: \$(date)"
    echo ""

    # Verify BiG-SLiCE installation
    echo "BiG-SLiCE version:"
    bigslice --version

    # Download models via BiG-SLiCE's built-in script
    echo ""
    echo "Downloading HMM models..."
    download_bigslice_hmmdb

    # Copy models from conda environment to output directory
    MODELS_SOURCE="\${CONDA_PREFIX}/bin/bigslice-models"

    echo ""
    echo "Verifying download..."
    if [ ! -d "\${MODELS_SOURCE}" ]; then
        echo "ERROR: BiG-SLiCE models not found"
        echo "Expected: \${MODELS_SOURCE}"
        exit 1
    fi

    # Create output directory and copy models
    mkdir -p hmmdb
    cp -r "\${MODELS_SOURCE}"/* hmmdb/

    # Verify copied files
    if [ ! -d "hmmdb/biosynthetic_pfams" ]; then
        echo "ERROR: Models not copied correctly"
        exit 1
    fi

    echo ""
    echo "Models directory contents:"
    ls -la hmmdb/
    echo ""
    echo "Biosynthetic Pfams (sample):"
    ls -la hmmdb/biosynthetic_pfams/ | head -10

    echo ""
    echo "=============================================="
    echo "BiG-SLiCE models ready"
    echo "Completed at: \$(date)"
    echo "=============================================="
    """
}
