process GTDBTK_CLASSIFY {
    tag "${taxon}"
    label 'process_high_memory'
    publishDir "${params.outdir}/gtdbtk_results/${Utils.sanitizeTaxon(taxon)}", mode: 'copy'
    cache 'lenient'  // GTDB-Tk is memory-intensive - use lenient caching

    input:
    val taxon
    path fasta_files
    path gtdbtk_db

    output:
    path "gtdbtk_output", emit: output_dir
    // GTDB-Tk v2.x outputs - tree files include class index (e.g., .1.tree)
    path "gtdbtk_output/classify/gtdbtk.bac120.classify.tree.*.tree", emit: bacterial_tree, optional: true
    path "gtdbtk_output/classify/gtdbtk.ar53.classify.tree.*.tree", emit: archaeal_tree, optional: true
    // Summary files are at top-level in GTDB-Tk v2.x output
    path "gtdbtk_output/gtdbtk.bac120.summary.tsv", emit: bacterial_summary, optional: true
    path "gtdbtk_output/gtdbtk.ar53.summary.tsv", emit: archaeal_summary, optional: true

    script:
    def cpus = params.gtdbtk_cpus ?: task.cpus ?: 8  // Use dedicated param, fallback to task.cpus
    def pplacer_cpus = params.gtdbtk_pplacer_cpus ?: 1  // pplacer is memory-heavy, use 1 by default
    def skip_ani = "--skip_ani_screen"
    def min_perc_aa = params.gtdbtk_min_perc_aa ?: 10
    """
    echo "=============================================="
    echo "GTDB-Tk Classification"
    echo "=============================================="
    echo "Taxon: ${taxon}"
    echo "CPUs: ${cpus}"
    echo "Pplacer CPUs: ${pplacer_cpus}"
    echo ""

    # Set GTDB-Tk database path
    RELEASE_DIR=\$(ls -d ${gtdbtk_db}/release* 2>/dev/null | head -1)
    if [ -z "\$RELEASE_DIR" ]; then
        export GTDBTK_DATA_PATH="${gtdbtk_db}"
    else
        export GTDBTK_DATA_PATH="\$RELEASE_DIR"
    fi
    echo "Using GTDB-Tk database: \$GTDBTK_DATA_PATH"

    # Create input genome directory
    # Add unique prefix to avoid conflicts with GTDB reference genome names
    # This prevents "duplicate taxon labels" errors in tree generation
    mkdir -p genomes
    for f in ${fasta_files}; do
        BASENAME=\$(basename "\$f")
        cp "\$f" "genomes/usr_\${BASENAME}"
    done

    GENOME_COUNT=\$(ls genomes/*.fna 2>/dev/null | wc -l)
    echo "Processing \$GENOME_COUNT genomes"
    echo ""

    if [ "\$GENOME_COUNT" -eq 0 ]; then
        echo "ERROR: No FASTA files found"
        exit 1
    fi

    # Run GTDB-Tk classify workflow
    echo "Running GTDB-Tk classify_wf..."
    echo "This may take several hours for large datasets."
    echo ""

    # Create scratch directory to reduce peak memory usage
    mkdir -p scratch_tmp

    # Run GTDB-Tk - capture exit code but don't fail immediately
    # Tree generation can fail with duplicate taxon labels (known GTDB-Tk issue)
    # when user genomes match reference strain names exactly
    set +e
    gtdbtk classify_wf \\
        --genome_dir genomes \\
        --out_dir gtdbtk_output \\
        --extension fna \\
        --cpus ${cpus} \\
        --pplacer_cpus ${pplacer_cpus} \\
        --min_perc_aa ${min_perc_aa} \\
        --scratch_dir scratch_tmp \\
        ${skip_ani} \\
        --prefix gtdbtk
    GTDBTK_EXIT=\$?
    set -e

    echo ""
    echo "=============================================="
    echo "GTDB-Tk classification complete (exit code: \$GTDBTK_EXIT)"
    echo "=============================================="

    # Check if we have summary files (main output) even if tree generation failed
    if [ \$GTDBTK_EXIT -ne 0 ]; then
        if [ -f gtdbtk_output/gtdbtk.bac120.summary.tsv ] || [ -f gtdbtk_output/gtdbtk.ar53.summary.tsv ]; then
            echo "WARNING: GTDB-Tk had errors but summary files were generated."
            echo "This often occurs due to duplicate taxon labels in tree generation."
            echo "Classification results are still valid."
        else
            echo "ERROR: GTDB-Tk failed and no summary files were generated."
            exit 1
        fi
    fi

    # Check outputs
    if [ -f gtdbtk_output/gtdbtk.bac120.summary.tsv ]; then
        BAC_COUNT=\$(tail -n +2 gtdbtk_output/gtdbtk.bac120.summary.tsv | wc -l)
        echo "Bacterial genomes classified: \$BAC_COUNT"
    fi
    if [ -f gtdbtk_output/gtdbtk.ar53.summary.tsv ]; then
        AR_COUNT=\$(tail -n +2 gtdbtk_output/gtdbtk.ar53.summary.tsv | wc -l)
        echo "Archaeal genomes classified: \$AR_COUNT"
    fi

    # List tree files
    echo ""
    echo "Phylogenetic tree files:"
    ls -la gtdbtk_output/classify/*.tree 2>/dev/null || echo "  No tree files found"
    """
}
