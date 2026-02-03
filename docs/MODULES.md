# ClusterQuest Module Reference

This document provides a detailed summary of each module in the ClusterQuest pipeline.

## Pipeline Architecture

```
main.nf
├── DOWNLOAD_GENOMES (subworkflow)
│   ├── NCBI_DATASETS_DOWNLOAD
│   ├── CREATE_NAME_MAP
│   ├── RENAME_GENOMES
│   └── EXTRACT_TAXONOMY
│
└── BGC_ANALYSIS (subworkflow)
    ├── ANTISMASH_ANALYSIS (subworkflow)
    │   ├── GET_ANTISMASH_VERSION
    │   ├── CHECK_ANTISMASH_REUSE
    │   ├── ANTISMASH
    │   └── COPY_ANTISMASH_RESULT
    │
    ├── Region Analysis
    │   ├── COUNT_REGIONS
    │   ├── TABULATE_REGIONS
    │   └── AGGREGATE_TAXONOMY
    │
    ├── CLUSTERING (subworkflow)
    │   ├── BIGSCAPE
    │   ├── BIGSLICE
    │   ├── EXTRACT_CLUSTERING_STATS
    │   └── EXTRACT_GCF_REPRESENTATIVES
    │
    ├── PHYLOGENY (subworkflow)
    │   ├── CHECK_GTDBTK_REUSE
    │   ├── GENBANK_TO_FASTA
    │   ├── GTDBTK_CLASSIFY
    │   └── FILTER_GTDBTK_RESULTS
    │
    ├── VISUALIZE_RESULTS
    └── COLLECT_VERSIONS
```

---

## Database Modules

### DOWNLOAD_ANTISMASH_DBS
**Location:** `modules/databases/download_antismash_dbs.nf`

Downloads antiSMASH databases required for BGC detection and annotation.

| Property | Value |
|----------|-------|
| Label | `process_low`, `retry_on_error` |
| Output | `db_dir` - antiSMASH database directory (~50 GB) |
| Cache | `storeDir` - downloaded once and reused |

### DOWNLOAD_PFAM
**Location:** `modules/databases/download_pfam.nf`

Downloads Pfam HMM profiles used by BiG-SCAPE for domain annotation.

| Property | Value |
|----------|-------|
| Label | `process_low`, `retry_on_error` |
| Output | `pfam_db` - Pfam-A.hmm database |
| Cache | `storeDir` |

### DOWNLOAD_GTDBTK_DB
**Location:** `modules/databases/download_gtdbtk_db.nf`

Downloads GTDB-Tk reference database for phylogenetic placement.

| Property | Value |
|----------|-------|
| Label | `process_low`, `retry_on_error` |
| Output | `db_dir` - GTDB-Tk database (~140 GB) |
| Cache | `storeDir` |

### DOWNLOAD_TAXONKIT_DB
**Location:** `modules/databases/download_taxonkit_db.nf`

Downloads NCBI TaxDump for taxonomy processing.

| Property | Value |
|----------|-------|
| Label | `process_low`, `retry_on_error` |
| Output | `taxdump_dir` - NCBI taxonomy database |
| Cache | `storeDir` |

### DOWNLOAD_BIGSLICE_MODELS
**Location:** `modules/databases/download_bigslice_models.nf`

Downloads BiG-SLiCE machine learning models.

| Property | Value |
|----------|-------|
| Label | `process_low`, `retry_on_error` |
| Output | `models_dir` - BiG-SLiCE models |
| Cache | `storeDir` |

---

## Genome Processing Modules

### NCBI_DATASETS_DOWNLOAD
**Location:** `modules/genome/ncbi_datasets_download.nf`

Downloads genomes from NCBI using the datasets CLI for a given taxon.

| Property | Value |
|----------|-------|
| Label | `process_low`, `retry_on_error` |
| Input | `taxon` - NCBI taxon name |
| Output | `genomes` - GenBank files, `assembly_info` - metadata |

**Features:**
- Uses `--dehydrated` download for efficient retrieval
- Validates downloads with checksum verification
- Handles large taxon sets (1000+ genomes)

### CREATE_NAME_MAP
**Location:** `modules/genome/create_name_map.nf`

Creates a mapping from NCBI assembly IDs to standardized genome names.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | `assembly_info` - NCBI metadata |
| Output | `name_map.json` - assembly ID → genome name |

**Naming rules:**
- Format: `Genus_species_strain` (e.g., `Streptomyces_coelicolor_A32`)
- Handles case-insensitive duplicates (appends numeric suffix)
- Sanitizes special characters for filesystem compatibility

### RENAME_GENOMES
**Location:** `modules/genome/rename_genomes_parallel.nf`

Renames genome files using the standardized naming scheme.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | `genome_pairs`, `name_map` |
| Output | `renamed_genome` - standardized .gbff files |

### GENBANK_TO_FASTA
**Location:** `modules/genome/genbank_to_fasta.nf`

Converts GenBank files to FASTA format for GTDB-Tk.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | GenBank file |
| Output | FASTA file (.fna) |

---

## BGC Analysis Modules

### GET_ANTISMASH_VERSION
**Location:** `modules/analysis/antismash.nf`

Retrieves the installed antiSMASH version for result tracking.

| Property | Value |
|----------|-------|
| Label | `process_local` |
| Output | Version string (e.g., "7.1.0") |

### ANTISMASH
**Location:** `modules/analysis/antismash.nf`

Runs antiSMASH BGC detection on a genome.

| Property | Value |
|----------|-------|
| Label | `process_medium`, `tolerant` |
| Input | `genome` - GenBank file, `db_dir` - databases |
| Output | `result_dir` - antiSMASH output directory |

**Features:**
- Configurable analyses: KnownClusterBlast, ClusterBlast, ClusterCompare
- Always enables clusterhmmer and tigrfam for domain analysis
- Writes `.antismash_meta` file for version/params tracking
- Skips genomes that already have results in publishDir

### CHECK_ANTISMASH_REUSE
**Location:** `modules/analysis/check_antismash_reuse.nf`

Checks if antiSMASH results can be reused from a previous taxon run.

| Property | Value |
|----------|-------|
| Label | `process_local` |
| Input | `genome`, `reuse_taxon`, `version`, `params_hash` |
| Output | `check_result` tuple: (genome, status, path) |

**Reuse criteria:**
- Result directory exists with valid .json output
- antiSMASH version matches (major.minor)
- Parameter hash matches current configuration

### COPY_ANTISMASH_RESULT
**Location:** `modules/analysis/check_antismash_reuse.nf`

Copies reused antiSMASH results to current taxon directory.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | Existing result path |
| Output | Copied result directory |

### COUNT_REGIONS
**Location:** `modules/analysis/count_regions.nf`

Counts BGC regions per genome from antiSMASH results.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | antiSMASH result directories |
| Output | `region_counts.tsv` |

**Options:**
- `--by_contig`: Count per contig instead of per genome
- `--split_hybrids`: Split hybrid BGC types

### TABULATE_REGIONS
**Location:** `modules/analysis/tabulate_regions.nf`

Creates detailed tabulation of all BGC regions.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | antiSMASH result directories |
| Output | `region_tabulation.tsv` |

**Columns include:**
- Genome, contig, region, product type
- Location (start, end, strand)
- KnownClusterBlast hits and similarity scores

### EXTRACT_TAXONOMY
**Location:** `modules/analysis/extract_taxonomy.nf`

Extracts taxonomy information from NCBI metadata.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | Assembly report, TaxDump database |
| Output | `taxonomy_map.json` |

### AGGREGATE_TAXONOMY
**Location:** `modules/analysis/aggregate_taxonomy.nf`

Aggregates BGC counts by taxonomic level.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | `taxonomy_map`, `region_counts` |
| Output | `taxonomy_tree.json` |

---

## Clustering Modules

### BIGSCAPE
**Location:** `modules/clustering/bigscape.nf`

Runs BiG-SCAPE network-based BGC clustering.

| Property | Value |
|----------|-------|
| Label | `process_high`, `tolerant` |
| Input | antiSMASH results, Pfam database |
| Output | `bigscape_dir`, `bigscape_db` (SQLite) |

**Features:**
- Configurable distance cutoffs (default: 0.30)
- Optional MIBiG integration
- Alignment modes: auto, global, glocal
- Generates SQLite database for downstream analysis

### BIGSLICE
**Location:** `modules/clustering/bigslice.nf`

Runs BiG-SLiCE ML-based BGC clustering.

| Property | Value |
|----------|-------|
| Label | `process_high` |
| Input | antiSMASH results, BiG-SLiCE models |
| Output | `bigslice_dir` |

**Features:**
- Configurable distance threshold
- Supports query mode against existing databases
- Integrates taxonomy information

### EXTRACT_CLUSTERING_STATS
**Location:** `modules/clustering/extract_clustering_stats.nf`

Extracts statistics from BiG-SCAPE or BiG-SLiCE results.

| Property | Value |
|----------|-------|
| Label | `process_low` |
| Input | Clustering output directory, tool name |
| Output | `{tool}_statistics.json` |

**Statistics include:**
- Total BGCs, GCFs, singletons
- BGC type distribution
- Clustering coverage

### EXTRACT_GCF_REPRESENTATIVES
**Location:** `modules/clustering/extract_gcf_representatives.nf`

Extracts representative BGCs for each GCF with gene diagrams.

| Property | Value |
|----------|-------|
| Label | `process_medium` |
| Input | BiG-SCAPE results, antiSMASH results |
| Output | `gcf_representatives.json` |

**Output includes:**
- Representative BGC per GCF
- Gene arrow diagrams with functional annotations
- KnownClusterBlast hit information

---

## Phylogeny Modules

### CHECK_GTDBTK_REUSE
**Location:** `modules/phylogeny/check_gtdbtk_reuse.nf`

Checks if GTDB-Tk results can be reused from a previous run.

| Property | Value |
|----------|-------|
| Label | `process_local` |
| Input | Genome list, reuse taxon |
| Output | Status (REUSE/RUN), paths to reuse files |

**Reuse criteria:**
- ALL current genomes must exist in previous results
- Summary TSV and tree file must be present

### GTDBTK_CLASSIFY
**Location:** `modules/phylogeny/gtdbtk.nf`

Runs GTDB-Tk classify workflow for phylogenetic placement.

| Property | Value |
|----------|-------|
| Label | `process_high_memory` |
| Input | FASTA files, GTDB-Tk database |
| Output | `bacterial_summary`, `bacterial_tree` |

**Resource requirements:**
- 56-64 GB RAM (pplacer step)
- ~140 GB disk for database

### FILTER_GTDBTK_RESULTS
**Location:** `modules/phylogeny/check_gtdbtk_reuse.nf`

Filters and prunes GTDB-Tk results for a genome subset.

| Property | Value |
|----------|-------|
| Label | `process_medium` |
| Input | Genome list, source summary, source tree |
| Output | Filtered summary TSV, pruned Newick tree |

---

## Visualization Modules

### VISUALIZE_RESULTS
**Location:** `modules/visualization/visualize_results.nf`

Generates the interactive HTML report.

| Property | Value |
|----------|-------|
| Label | `process_medium` |
| Input | All analysis outputs |
| Output | `bgc_report.html`, individual genome pages |

**Report sections:**
- Overview: BGC statistics, donut chart, rarefaction curve
- Taxonomy: Interactive taxonomy tree with BGC counts
- BGC Distribution: GCF × taxonomy heatmap
- Genomes: Searchable genome table
- Clustering: GCF visualization with gene diagrams
- Novel BGCs: Potentially novel clusters
- KCB Hits: Known cluster matches
- Resources: Pipeline timing and resource usage

### COLLECT_VERSIONS
**Location:** `modules/utilities/collect_versions.nf`

Collects software versions from conda environments.

| Property | Value |
|----------|-------|
| Label | `process_local` |
| Output | `software_versions.json` |

---

## Configuration Files

### conf/conda.config
Centralized conda environment specifications for all processes.

### conf/labels.config
Process labels for resource allocation:

| Label | CPUs | Memory | Time |
|-------|------|--------|------|
| `process_local` | 1 | 1 GB | 1h |
| `process_low` | 1 | 2 GB | 1h |
| `process_medium` | 4 | 8 GB | 4h |
| `process_high` | 8 | 32 GB | 8h |
| `process_high_memory` | 8 | 128 GB | 24h |

Error handling labels:
- `tolerant`: `errorStrategy = 'ignore'`
- `retry_on_error`: Retry on transient failures (up to 3 times)

---

## Utility Library

### lib/Utils.groovy

Helper functions used across modules:

| Function | Description |
|----------|-------------|
| `sanitizeTaxon(name)` | Sanitize taxon for filesystem paths |
| `antismashParamsHash(params)` | Generate MD5 hash of antiSMASH parameters |
| `buildReusePath(params, projectDir, tool, taxon, subPath)` | Build absolute path for result reuse |
| `isValidInput(input)` | Check if input is valid (not placeholder) |
