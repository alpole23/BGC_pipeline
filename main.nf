#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

/**
 * Check if a clustering method is enabled.
 */
def clusteringEnabled(method) {
    params.clustering == method || params.clustering == "both"
}

/**
 * Create placeholder channel for optional inputs.
 */
def placeholder(name) {
    Channel.value(file(name))
}

/**
 * Validate pipeline parameters.
 */
def validateParams() {
    def errors = []

    // Validate workflow
    def validWorkflows = ['download', 'bgc_analysis', 'full']
    if (!(params.workflow in validWorkflows)) {
        errors << "Invalid workflow '${params.workflow}'. Valid options: ${validWorkflows.join(', ')}"
    }

    // Validate clustering
    def validClustering = ['none', 'bigscape', 'bigslice', 'both']
    if (!(params.clustering in validClustering)) {
        errors << "Invalid clustering '${params.clustering}'. Valid options: ${validClustering.join(', ')}"
    }

    // Validate bigscape alignment mode
    def validAlignmentModes = ['auto', 'global', 'glocal']
    if (!(params.bigscape_alignment_mode in validAlignmentModes)) {
        errors << "Invalid bigscape_alignment_mode '${params.bigscape_alignment_mode}'. Valid options: ${validAlignmentModes.join(', ')}"
    }

    // Validate bgc_analysis workflow has required input
    if (params.workflow == 'bgc_analysis' && (!params.input_genomes || params.input_genomes == 'null')) {
        errors << "params.input_genomes must be specified for 'bgc_analysis' workflow"
    }

    // Validate reuse_antismash_from if specified
    if (params.reuse_antismash_from) {
        def reuse_dir = "${params.outdir}/antismash_results/${Utils.sanitizeTaxon(params.reuse_antismash_from)}"
        if (!file(reuse_dir).exists()) {
            errors << "reuse_antismash_from directory does not exist: ${reuse_dir}"
        }
    }

    // Validate reuse_gtdbtk_from if specified
    if (params.reuse_gtdbtk_from) {
        def reuse_dir = "${params.outdir}/gtdbtk_results/${Utils.sanitizeTaxon(params.reuse_gtdbtk_from)}/gtdbtk_output"
        if (!file(reuse_dir).exists()) {
            errors << "reuse_gtdbtk_from directory does not exist: ${reuse_dir}"
        }
    }

    // Report errors
    if (errors) {
        log.error "=" * 60
        log.error "PARAMETER VALIDATION FAILED"
        log.error "=" * 60
        errors.each { log.error "  - ${it}" }
        log.error "=" * 60
        error "Please fix the above parameter errors and try again."
    }

    // Warnings
    if (params.run_gtdbtk) {
        log.warn "GTDB-Tk is enabled. This requires ~140 GB disk space and ~56-64 GB RAM."
    }

    if (params.reuse_antismash_from) {
        log.info "antiSMASH result reuse enabled from taxon: ${params.reuse_antismash_from}"
    }

    if (params.reuse_gtdbtk_from) {
        log.info "GTDB-Tk result reuse enabled from taxon: ${params.reuse_gtdbtk_from}"
    }
}

// Run validation
validateParams()

// =============================================================================
// MODULE IMPORTS
// =============================================================================

// Database download modules
include { DOWNLOAD_PFAM } from './modules/databases/download_pfam'
include { DOWNLOAD_ANTISMASH_DBS } from './modules/databases/download_antismash_dbs'
include { DOWNLOAD_TAXONKIT_DB } from './modules/databases/download_taxonkit_db'
include { DOWNLOAD_GTDBTK_DB } from './modules/databases/download_gtdbtk_db'
include { DOWNLOAD_BIGSLICE_MODELS } from './modules/databases/download_bigslice_models'

// Genome processing modules
include { NCBI_DATASETS_DOWNLOAD } from './modules/genome/ncbi_datasets_download'
include { CREATE_NAME_MAP } from './modules/genome/create_name_map'
include { RENAME_GENOMES } from './modules/genome/rename_genomes_parallel'
include { GENBANK_TO_FASTA } from './modules/genome/genbank_to_fasta'

// BGC analysis modules
include { GET_ANTISMASH_VERSION; ANTISMASH } from './modules/analysis/antismash'
include { CHECK_ANTISMASH_REUSE; COPY_ANTISMASH_RESULT } from './modules/analysis/check_antismash_reuse'
include { COUNT_REGIONS } from './modules/analysis/count_regions'
include { TABULATE_REGIONS } from './modules/analysis/tabulate_regions'
include { EXTRACT_TAXONOMY } from './modules/analysis/extract_taxonomy'
include { AGGREGATE_TAXONOMY } from './modules/analysis/aggregate_taxonomy'

// Clustering modules
include { BIGSCAPE } from './modules/clustering/bigscape'
include { BIGSLICE } from './modules/clustering/bigslice'
include { EXTRACT_CLUSTERING_STATS as EXTRACT_BIGSCAPE_STATS } from './modules/clustering/extract_clustering_stats'
include { EXTRACT_CLUSTERING_STATS as EXTRACT_BIGSLICE_STATS } from './modules/clustering/extract_clustering_stats'
include { EXTRACT_GCF_REPRESENTATIVES } from './modules/clustering/extract_gcf_representatives'

// Phylogenetic analysis modules
include { GTDBTK_CLASSIFY } from './modules/phylogeny/gtdbtk'
include { CHECK_GTDBTK_REUSE; FILTER_GTDBTK_RESULTS } from './modules/phylogeny/check_gtdbtk_reuse'

// Visualization modules
include { VISUALIZE_RESULTS } from './modules/visualization/visualize_results'

// Utility modules
include { COLLECT_VERSIONS } from './modules/utilities/collect_versions'

// =============================================================================
// SUBWORKFLOWS
// =============================================================================

/*
 * Subworkflow: Download and prepare genomes from NCBI
 */
workflow DOWNLOAD_GENOMES {
    take:
        taxon

    main:
        NCBI_DATASETS_DOWNLOAD(taxon)
        CREATE_NAME_MAP(taxon, NCBI_DATASETS_DOWNLOAD.out.assembly_info)
        DOWNLOAD_TAXONKIT_DB()

        // Prepare genome pairs (assembly_id, genome_file)
        genome_pairs = NCBI_DATASETS_DOWNLOAD.out.genomes
            .flatten()
            .map { gbff -> tuple(gbff.parent.name, gbff) }

        RENAME_GENOMES(taxon, genome_pairs, CREATE_NAME_MAP.out.name_map)

        EXTRACT_TAXONOMY(
            taxon,
            NCBI_DATASETS_DOWNLOAD.out.assembly_data_report,
            NCBI_DATASETS_DOWNLOAD.out.taxonomy_report,
            DOWNLOAD_TAXONKIT_DB.out.taxdump_dir
        )

    emit:
        renamed_genomes      = RENAME_GENOMES.out.renamed_genome
        assembly_info        = NCBI_DATASETS_DOWNLOAD.out.assembly_info
        name_map             = CREATE_NAME_MAP.out.name_map
        taxonomy_map         = EXTRACT_TAXONOMY.out.taxonomy_map
        assembly_data_report = NCBI_DATASETS_DOWNLOAD.out.assembly_data_report
        taxonomy_report      = NCBI_DATASETS_DOWNLOAD.out.taxonomy_report
}

/*
 * Subworkflow: Run antiSMASH on genomes (with optional result reuse)
 */
workflow ANTISMASH_ANALYSIS {
    take:
        taxon
        renamed_genomes

    main:
        DOWNLOAD_ANTISMASH_DBS()

        // Get antiSMASH version for tracking
        GET_ANTISMASH_VERSION()
        antismash_version = GET_ANTISMASH_VERSION.out.version

        // Generate hash of current antiSMASH parameters for tracking
        antismash_params_hash = Utils.antismashParamsHash(params)

        if (params.reuse_antismash_from) {
            // === REUSE MODE ===
            CHECK_ANTISMASH_REUSE(
                taxon,
                params.reuse_antismash_from,
                renamed_genomes,
                antismash_version,
                antismash_params_hash
            )

            // Split: genomes to run vs genomes to reuse
            genomes_to_run = CHECK_ANTISMASH_REUSE.out.check_result
                .filter { genome, status, path -> status == "RUN" }
                .map { genome, status, path -> genome }

            genomes_to_reuse = CHECK_ANTISMASH_REUSE.out.check_result
                .filter { genome, status, path -> status == "REUSE" }
                .map { genome, status, path -> tuple(genome.baseName, file(path)) }

            // Run antiSMASH on genomes that need it
            ANTISMASH(taxon, genomes_to_run, DOWNLOAD_ANTISMASH_DBS.out.db_dir, antismash_version, antismash_params_hash)

            // Copy reused results
            COPY_ANTISMASH_RESULT(
                taxon,
                genomes_to_reuse.map { it[0] },
                genomes_to_reuse.map { it[1] }
            )

            // Combine all results
            antismash_results = ANTISMASH.out.result_dir
                .mix(COPY_ANTISMASH_RESULT.out.result_dir)
                .collect()
        } else {
            // === NORMAL MODE ===
            ANTISMASH(taxon, renamed_genomes, DOWNLOAD_ANTISMASH_DBS.out.db_dir, antismash_version, antismash_params_hash)
            antismash_results = ANTISMASH.out.result_dir.collect()
        }

    emit:
        results = antismash_results
}

/*
 * Subworkflow: BiG-SCAPE and BiG-SLiCE clustering
 */
workflow CLUSTERING {
    take:
        taxon
        antismash_results
        taxonomy_map
        name_map
        tabulation

    main:
        bigscape_stats_ch = placeholder('NO_BIGSCAPE_STATS')
        bigscape_db_ch = placeholder('NO_BIGSCAPE_DB')
        bigscape_dir_ch = placeholder('NO_BIGSCAPE_DIR')
        bigslice_stats_ch = placeholder('NO_BIGSLICE_STATS')
        gcf_data_ch = placeholder('NO_GCF_DATA')

        if (clusteringEnabled("bigscape")) {
            DOWNLOAD_PFAM()
            BIGSCAPE(taxon, antismash_results, DOWNLOAD_PFAM.out.pfam_db)
            EXTRACT_BIGSCAPE_STATS(taxon, "bigscape", BIGSCAPE.out.bigscape_dir)
            bigscape_stats_ch = EXTRACT_BIGSCAPE_STATS.out.stats_json
            bigscape_db_ch = BIGSCAPE.out.bigscape_db
            bigscape_dir_ch = BIGSCAPE.out.bigscape_dir

            // Extract GCF representatives (needs tabulation for KCB hit lookup)
            if (tabulation.name != 'NO_TABULATION') {
                EXTRACT_GCF_REPRESENTATIVES(taxon, bigscape_dir_ch, antismash_results, tabulation)
                gcf_data_ch = EXTRACT_GCF_REPRESENTATIVES.out.gcf_data
            }
        }

        if (clusteringEnabled("bigslice")) {
            DOWNLOAD_BIGSLICE_MODELS()
            BIGSLICE(taxon, antismash_results, DOWNLOAD_BIGSLICE_MODELS.out.models_dir, taxonomy_map, name_map)
            EXTRACT_BIGSLICE_STATS(taxon, "bigslice", BIGSLICE.out.bigslice_dir)
            bigslice_stats_ch = EXTRACT_BIGSLICE_STATS.out.stats_json
        }

    emit:
        bigscape_stats = bigscape_stats_ch
        bigscape_db    = bigscape_db_ch
        bigslice_stats = bigslice_stats_ch
        gcf_data       = gcf_data_ch
}

/*
 * Subworkflow: GTDB-Tk phylogenetic classification (with optional result reuse)
 */
workflow PHYLOGENY {
    take:
        taxon
        renamed_genomes
        counts

    main:
        phylo_tree_ch = placeholder('NO_PHYLO_TREE')
        gtdbtk_summary_ch = placeholder('NO_GTDBTK_SUMMARY')

        if (params.run_gtdbtk) {
            // Determine which genomes to process
            if (params.gtdbtk_bgc_genomes_only) {
                // Filter to genomes with BGCs
                genomes_with_bgcs_ch = counts
                    .splitCsv(header: true, sep: '\t', skip: 1)
                    .filter { row -> (row.total_count ?: '0').toInteger() > 0 }
                    .map { row -> tuple(row.record, true) }

                renamed_genomes_tuples = renamed_genomes
                    .map { genome -> tuple(genome.name, genome) }

                genomes_for_gtdbtk = renamed_genomes_tuples
                    .join(genomes_with_bgcs_ch)
                    .map { name, genome, flag -> genome }
            } else {
                genomes_for_gtdbtk = renamed_genomes
            }

            // Convert GenBank to FASTA
            GENBANK_TO_FASTA(genomes_for_gtdbtk)
            fasta_files = GENBANK_TO_FASTA.out.fasta.collect()

            if (params.reuse_gtdbtk_from) {
                // === GTDB-Tk REUSE MODE ===
                genome_list_ch = GENBANK_TO_FASTA.out.fasta
                    .map { it.toString() }
                    .collectFile(name: 'genome_list.txt', newLine: true)

                CHECK_GTDBTK_REUSE(taxon, params.reuse_gtdbtk_from, genome_list_ch)

                check_result = CHECK_GTDBTK_REUSE.out.check_result
                    .branch {
                        reuse: it[0] == "REUSE"
                        run: it[0] == "RUN"
                    }

                // REUSE path
                FILTER_GTDBTK_RESULTS(
                    taxon,
                    genome_list_ch,
                    check_result.reuse.map { it[1] },
                    check_result.reuse.map { it[2] }
                )

                // RUN path
                DOWNLOAD_GTDBTK_DB()
                fasta_for_fresh_run = check_result.run
                    .combine(fasta_files)
                    .map { status, summary, tree, files -> files }
                GTDBTK_CLASSIFY(taxon, fasta_for_fresh_run, DOWNLOAD_GTDBTK_DB.out.db_dir)

                // Combine outputs
                phylo_tree_ch = FILTER_GTDBTK_RESULTS.out.bacterial_tree
                    .mix(GTDBTK_CLASSIFY.out.bacterial_tree)
                    .ifEmpty(file('NO_PHYLO_TREE'))
                gtdbtk_summary_ch = FILTER_GTDBTK_RESULTS.out.bacterial_summary
                    .mix(GTDBTK_CLASSIFY.out.bacterial_summary)
                    .ifEmpty(file('NO_GTDBTK_SUMMARY'))
            } else {
                // === GTDB-Tk NORMAL MODE ===
                DOWNLOAD_GTDBTK_DB()
                GTDBTK_CLASSIFY(taxon, fasta_files, DOWNLOAD_GTDBTK_DB.out.db_dir)

                phylo_tree_ch = GTDBTK_CLASSIFY.out.bacterial_tree.ifEmpty(file('NO_PHYLO_TREE'))
                gtdbtk_summary_ch = GTDBTK_CLASSIFY.out.bacterial_summary.ifEmpty(file('NO_GTDBTK_SUMMARY'))
            }
        }

    emit:
        tree    = phylo_tree_ch
        summary = gtdbtk_summary_ch
}

// =============================================================================
// MAIN WORKFLOWS
// =============================================================================

/*
 * Workflow: Complete BGC detection and analysis pipeline
 */
workflow BGC_ANALYSIS {
    take:
        taxon
        renamed_genomes
        assembly_info
        name_map
        taxonomy_map

    main:
        // --- BGC Detection ---
        ANTISMASH_ANALYSIS(taxon, renamed_genomes)
        antismash_results = ANTISMASH_ANALYSIS.out.results

        // --- Region Analysis ---
        counts_ch = placeholder('NO_COUNTS')
        tabulation_ch = placeholder('NO_TABULATION')
        taxonomy_tree_ch = placeholder('NO_TAXONOMY_TREE')

        if (params.run_analysis) {
            COUNT_REGIONS(taxon, antismash_results)
            counts_ch = COUNT_REGIONS.out.counts

            AGGREGATE_TAXONOMY(taxon, taxonomy_map, COUNT_REGIONS.out.counts, name_map)
            taxonomy_tree_ch = AGGREGATE_TAXONOMY.out.taxonomy_tree

            TABULATE_REGIONS(taxon, antismash_results)
            tabulation_ch = TABULATE_REGIONS.out.tabulation
        }

        // --- Clustering ---
        CLUSTERING(taxon, antismash_results, taxonomy_map, name_map, tabulation_ch)

        // --- Phylogenetic Analysis ---
        PHYLOGENY(taxon, renamed_genomes, counts_ch)

        // --- Visualization ---
        if (params.run_analysis) {
            COLLECT_VERSIONS(antismash_results, CLUSTERING.out.bigscape_db, PHYLOGENY.out.summary)
            versions_ch = COLLECT_VERSIONS.out.versions

            trace_file_ch = file("${params.outdir}/pipeline_info/pipeline_trace.tsv").exists()
                ? Channel.value(file("${params.outdir}/pipeline_info/pipeline_trace.tsv"))
                : Channel.value(file('NO_TRACE_FILE'))

            VISUALIZE_RESULTS(
                taxon,
                counts_ch,
                tabulation_ch,
                assembly_info,
                name_map,
                taxonomy_map,
                taxonomy_tree_ch,
                CLUSTERING.out.bigslice_stats,
                CLUSTERING.out.bigscape_stats,
                CLUSTERING.out.bigscape_db,
                CLUSTERING.out.gcf_data,
                PHYLOGENY.out.tree,
                PHYLOGENY.out.summary,
                trace_file_ch,
                versions_ch
            )
        }
}

/*
 * Main entry point
 */
workflow {
    def taxon_dir = Utils.sanitizeTaxon(params.taxon)

    switch (params.workflow) {
        case "download":
            DOWNLOAD_GENOMES(params.taxon)
            break

        case "bgc_analysis":
            // Load from previous download run
            def base_dir = "${params.outdir}/ncbi_genomes/${taxon_dir}"
            def results_dir = "${params.outdir}/main_analysis_results/${taxon_dir}"

            BGC_ANALYSIS(
                params.taxon,
                Channel.fromPath("${params.input_genomes}/*.gbff"),
                Channel.fromPath("${base_dir}/ncbi_dataset/data/assembly_info_table.txt"),
                Channel.fromPath("${base_dir}/name_map.json"),
                Channel.fromPath("${results_dir}/taxonomy_map.json")
            )
            break

        case "full":
            DOWNLOAD_GENOMES(params.taxon)
            BGC_ANALYSIS(
                params.taxon,
                DOWNLOAD_GENOMES.out.renamed_genomes,
                DOWNLOAD_GENOMES.out.assembly_info,
                DOWNLOAD_GENOMES.out.name_map,
                DOWNLOAD_GENOMES.out.taxonomy_map
            )
            break

        default:
            error "ERROR: Invalid workflow '${params.workflow}'. Valid options: 'download', 'bgc_analysis', 'full'"
    }
}
