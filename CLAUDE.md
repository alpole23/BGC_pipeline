# CLAUDE.md

Nextflow pipeline for analyzing phosphonate biosynthetic gene clusters (BGCs) in bacterial genomes using antiSMASH with optional BiG-SCAPE clustering and GTDB-Tk phylogenetic analysis.

## Environment Setup

```bash
conda activate nextflow    # Activate conda environment before running
```

## Quick Start

```bash
nextflow run main.nf --taxon "Pantoea ananatis"           # Full pipeline
nextflow run main.nf -resume                               # Resume previous run
nextflow run main.nf --workflow download --taxon "Streptomyces coelicolor"
nextflow run main.nf --clustering bigscape                 # With clustering
nextflow run main.nf -profile slurm                        # HPC execution
```

## Cross-Taxon Result Reuse

When analyzing a taxon that's a subset of a previously analyzed taxon, you can reuse existing antiSMASH and GTDB-Tk results to avoid redundant computation.

### Usage

```bash
# First run on broad taxon (e.g., family level)
nextflow run main.nf --taxon "Erwiniaceae"

# Later, run on subset taxon, reusing results
nextflow run main.nf --taxon "Pantoea" --reuse_antismash_from "Erwiniaceae" --reuse_gtdbtk_from "Erwiniaceae"
```

### antiSMASH Reuse

1. For each genome in the current run, the pipeline checks if results exist in the reuse directory
2. Results are reused if:
   - The antiSMASH version matches (major.minor, e.g., 7.1.x matches 7.1.y)
   - The parameter configuration matches (tracked via hash)
3. Genomes without existing results are processed normally
4. Each antiSMASH result includes a `.antismash_meta` file that stores version and params_hash

### GTDB-Tk Reuse

1. The pipeline checks if ALL genomes in the current run exist in the reuse results
2. If yes: The summary TSV is filtered and the phylogenetic tree is pruned to only include current genomes
3. If no: GTDB-Tk runs fresh on all current genomes (all-or-nothing approach)

This is more efficient than re-running GTDB-Tk, especially since the classify step (pplacer) is memory-intensive.

### Clustering

BiG-SCAPE always runs fresh for the current genome set, as clustering depends on the complete set of BGCs being analyzed together.

### When to Use

- Running on a genus after analyzing the family (e.g., Pantoea after Erwiniaceae)
- Re-running with different clustering parameters (antiSMASH results unchanged)
- Adding new genomes to a previous analysis

## Configuration

Parameters in `nextflow.config` are organized by subworkflow to make it easy to find relevant settings.

### Global

| Parameter | Default | Description |
|-----------|---------|-------------|
| `workflow` | "full" | Pipeline mode: `download`, `bgc_analysis`, or `full` |
| `outdir` | "results" | Output directory for all results |

### DOWNLOAD_GENOMES

Downloads and prepares bacterial genomes from NCBI.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `taxon` | "Pantoea ananatis" | NCBI taxon (species, genus, family, order, etc.) |

### ANTISMASH_ANALYSIS

BGC detection using antiSMASH.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_genomes` | null | Path to pre-downloaded genomes (for `bgc_analysis` workflow) |
| `reuse_antismash_from` | null | Taxon name to reuse antiSMASH results from |
| `antismash_minimal` | false | Minimal mode (faster, skips domain analysis) |
| `antismash_cb_general` | false | ClusterBlast: Compare vs antiSMASH DB |
| `antismash_cc_mibig` | false | ClusterCompare: Advanced MIBiG scoring |
| `antismash_smcog_trees` | false | Phylogenetic trees for BGC genes |

**Note:** Detection is hardcoded to phosphonate rule only (`--hmmdetection-limit-to-rule-names phosphonate`). `--cb-knownclusters`, `--clusterhmmer`, and `--tigrfam` are always enabled.

### Region Analysis

BGC counting, tabulation, and statistics.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_analysis` | true | Enable region analysis and visualization |
| `count_per_contig` | false | Count per contig (true) or per genome (false) |
| `split_hybrids` | false | Split hybrid types (T1PKS-NRPS → T1PKS + NRPS) |

### CLUSTERING

Gene Cluster Family (GCF) clustering.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `clustering` | "bigscape" | `none` or `bigscape` |

**BiG-SCAPE Options** (when `clustering = "bigscape"`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bigscape_cutoffs` | "0.30" | GCF distance threshold(s), comma-separated |
| `bigscape_alignment_mode` | "auto" | `auto`, `global`, or `glocal` |
| `bigscape_mibig_version` | "" | MIBiG version (e.g., "3.1") or "" to exclude |
| `bigscape_classify` | "category" | `""`, `category`, `class`, or `legacy` |
| `bigscape_include_singletons` | true | Include unclustered BGCs |
| `bigscape_mix` | false | Allow mixing BGC classes in same GCF |

### PHYLOGENY

Phylogenetic placement using GTDB-Tk. ⚠️ Requires ~140 GB disk and ~56-64 GB RAM.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_gtdbtk` | true | Enable phylogenetic analysis |
| `reuse_gtdbtk_from` | null | Taxon name to reuse GTDB-Tk results from |
| `gtdbtk_bgc_genomes_only` | true | Only analyze genomes with detected BGCs |
| `gtdbtk_cpus` | 8 | CPUs for identify/align steps |
| `gtdbtk_pplacer_cpus` | 1 | CPUs for pplacer (keep low, memory-bound) |
| `gtdbtk_min_perc_aa` | 10 | Minimum % amino acids in MSA |
| `gtdbtk_outgroup` | null | Outgroup pattern for tree rooting (e.g., "g__Escherichia") |

## Output Structure

```
results/
├── databases/                    # Cached databases (preserve these)
├── ncbi_genomes/${taxon}/
│   ├── ncbi_dataset/            # Raw NCBI download
│   ├── renamed_genomes/         # Standardized genome files
│   └── name_map.json            # Assembly ID to genome name mapping
├── antismash_results/${taxon}/  # Per-genome antiSMASH results
├── bigscape_results/${taxon}/   # BiG-SCAPE clustering output
│   ├── ${taxon}.db              # SQLite database with clustering results
│   └── gcf_representatives.json # GCF data with KCB hits and gene diagrams
├── gtdbtk_results/${taxon}/     # GTDB-Tk phylogenetic placement
├── pipeline_info/
│   ├── pipeline_trace.tsv       # Per-task timing and resource data
│   ├── pipeline_report.html     # Nextflow execution report
│   ├── pipeline_timeline.html   # Visual timeline
│   └── software_versions.json   # Tool versions
└── main_analysis_results/${taxon}/
    ├── region_counts.tsv        # BGC counts per genome
    ├── region_tabulation.tsv    # Detailed BGC information
    ├── taxonomy_map.json        # Genome taxonomy mapping
    ├── bgc_report.html          # Interactive HTML report
    ├── genomes/                 # Per-genome HTML pages
    ├── kcb_identification_chart.png
    ├── rarefaction_curve.png
    └── gcf_heatmap/             # GCF biosynthetic tree and heatmap outputs
        ├── gcf_biosynthetic_tree.png
        ├── gcf_biosynthetic_tree.svg
        ├── gcf_species_heatmap.png
        ├── gcf_species_heatmap.svg
        └── phosphonate_metadata.json
```

## SLURM HPC Execution

### Configuration

The pipeline uses process labels from `conf/labels.config` with SLURM-specific overrides in `nextflow.config`. The SLURM profile (`-profile slurm`) provides:

- Queue management: 200 concurrent jobs, 20 submissions/min
- High-memory queue for `process_high_memory` label (GTDB-Tk)
- BiG-SCAPE memory override: 64 GB (scales O(n²))
- Extended time limits for network-bound processes

Default resources are defined by labels and can be overridden in the SLURM profile.

### Running on SLURM

```bash
# Direct execution
nextflow run main.nf -profile slurm --taxon "Pantoea"

# Via sbatch (recommended for long runs)
sbatch submit_slurm.sh
```

## Benchmarking

### Purpose

Benchmark the pipeline to estimate runtime for large-scale analyses (e.g., 3 million genomes).

### Recommended Benchmark Size

| Sample Size | Use Case |
|-------------|----------|
| 100-200 | Good statistical power, captures variance |
| 300-500 | Excellent confidence intervals, tests SLURM scaling |
| 1000+ | Very accurate but excessive for benchmarking |

### Running a Benchmark

1. **Pre-download databases** (one-time cost, not included in benchmark):
   ```bash
   # Databases are cached in results/databases/
   # antiSMASH (~50GB), GTDB-Tk (~140GB), Pfam, TaxonKit
   ```

2. **Run benchmark with trace enabled** (already configured in nextflow.config):
   ```bash
   nextflow run main.nf -profile slurm --taxon "Pantoea"
   ```

3. **Analyze trace file**:
   ```bash
   # Trace file location: results/pipeline_info/pipeline_trace.tsv
   # Contains: task_id, name, status, realtime, cpus, memory, peak_rss, peak_vmem
   ```

### Extrapolation Formula

```
Total time = (mean_time_per_genome × total_genomes) / max_parallel_jobs

# With 95% confidence interval:
CI = mean ± (1.96 × std_dev / sqrt(n))
```

### Key Metrics

- **Wall time per genome**: Primary scaling factor (antiSMASH is bottleneck)
- **CPU hours per genome**: For HPC allocation requests
- **Peak memory per genome**: For SLURM memory allocation
- **Success/failure rate**: For planning retries

### Example: Pantoea Benchmark

- **Taxon**: Pantoea (~1500 genomes)
- **Configuration**: Full analysis (antiSMASH + BiG-SCAPE + GTDB-Tk)
- **Expected output**: Per-genome timing data for 3M genome extrapolation

## Module & Script Organization

```
main.nf                 # Main workflow (uses subworkflows)
nextflow.config         # Parameters and SLURM profile

conf/
├── conda.config        # Centralized conda environments by process
└── labels.config       # Process labels (resource allocations, error handling)

lib/
└── Utils.groovy        # Shared Groovy utilities (sanitizeTaxon, antismashParamsHash, buildReusePath)

modules/
├── databases/          # Database download processes (antiSMASH, GTDB-Tk, Pfam, etc.)
├── genome/             # Genome processing (NCBI download, rename, GenBank→FASTA)
├── analysis/           # BGC analysis (antiSMASH, counting, tabulation, reuse)
├── clustering/         # BiG-SCAPE clustering and stats extraction
├── phylogeny/          # GTDB-Tk classification (with reuse support)
├── visualization/      # HTML report generation
└── utilities/          # Version collection

scripts/
├── utils/              # Shared Python utilities
│   ├── constants.py      # BGC_COLORS, GENE_COLORS, KCB_THRESHOLDS
│   ├── parsers.py        # Duration, memory, timestamp parsing
│   └── antismash_parser.py  # antiSMASH JSON parsing
├── viz/                # Visualization modules
│   ├── charts.py         # KCB pie charts, BGC color utilities
│   ├── tree_viz.py       # Phylogenetic/taxonomy tree visualization
│   ├── tables.py         # Genome tables, statistics
│   ├── clustering.py     # BiG-SCAPE stats HTML
│   ├── taxonomy.py       # Interactive taxonomy tree
│   └── resources.py      # Resource usage visualization
├── taxonomy/           # Taxonomy processing scripts
├── genome/             # Genome processing scripts
├── clustering/         # Clustering statistics and GCF representative extraction
├── analysis/           # BGC counting and tabulation
├── phylogeny/          # GTDB-Tk result filtering
└── visualize_results.py  # Main visualization entry point
```

### Workflow Structure

The pipeline uses DSL2 subworkflows for modularity. Parameters in `nextflow.config` are organized to mirror this structure.

```
workflow (entry point)
│
├── DOWNLOAD_GENOMES          # Download and prepare genomes from NCBI
│   ├── NCBI_DATASETS_DOWNLOAD
│   ├── CREATE_NAME_MAP
│   ├── RENAME_GENOMES
│   └── EXTRACT_TAXONOMY
│
└── BGC_ANALYSIS              # Main analysis pipeline
    │
    ├── ANTISMASH_ANALYSIS    # BGC detection (with reuse support)
    │   ├── CHECK_ANTISMASH_REUSE
    │   ├── ANTISMASH
    │   └── COPY_ANTISMASH_RESULT
    │
    ├── Region Analysis       # BGC statistics
    │   ├── COUNT_REGIONS
    │   ├── TABULATE_REGIONS
    │   └── AGGREGATE_TAXONOMY
    │
    ├── CLUSTERING            # GCF clustering
    │   ├── BIGSCAPE
    │   ├── EXTRACT_CLUSTERING_STATS
    │   └── EXTRACT_GCF_REPRESENTATIVES
    │
    ├── PHYLOGENY             # GTDB-Tk (with reuse support)
    │   ├── CHECK_GTDBTK_REUSE
    │   ├── GTDBTK_CLASSIFY
    │   └── FILTER_GTDBTK_RESULTS
    │
    ├── GCF_BIOSYNTHETIC_TREE # GCF biosynthetic NJ tree (when bigscape enabled, runs before visualization)
    └── VISUALIZE_RESULTS     # HTML report generation (receives GCF tree PNG as input)
```

**Invoking subworkflows directly:**

```bash
# Run only download
nextflow run main.nf -entry DOWNLOAD_GENOMES --taxon "Pantoea"

# Run full pipeline (default)
nextflow run main.nf --taxon "Pantoea"
```

### Process Labels

Processes use labels for resource allocation and error handling:

| Label | CPUs | Memory | Description |
|-------|------|--------|-------------|
| `process_local` | 1 | 1 GB | Runs on head node |
| `process_low` | 1 | 2 GB | Light scripts |
| `process_medium` | 4 | 8 GB | antiSMASH, visualization |
| `process_high` | 8 | 32 GB | BiG-SCAPE |
| `process_high_memory` | 8 | 128 GB | GTDB-Tk pplacer |

Error handling labels:
- `tolerant`: Individual failures don't stop pipeline (per-genome processes)
- `retry_on_error`: Retry on transient errors (network downloads)

## Software Versions

Versions are dynamically collected from installed tools. Most use `--version` flag, but TaxonKit uses `version` subcommand.

| Tool | Conda Spec | Purpose |
|------|------------|---------|
| antiSMASH | bioconda::antismash | BGC detection |
| BiG-SCAPE | bioconda::bigscape | GCF clustering |
| GTDB-Tk | bioconda::gtdbtk | Phylogenetic placement |
| TaxonKit | bioconda::taxonkit | Taxonomy processing |

Version information is output to `results/pipeline_info/software_versions.json`.

## HTML Report Features

The interactive HTML report (`bgc_report.html`) includes:

### Tabs
The report uses 6 tabs:
- **Overview**: Summary statistics grid, rarefaction curve, pipeline resource usage (collapsible) and software versions
- **Phylogeny**: NCBI taxonomy tree + GTDB-Tk phylogenetic tree and BGC distribution
- **Genomes**: Searchable genome table with links to individual genome pages
- **GCF Analysis**: GCF biosynthetic NJ tree (embedded as base64), coupling enzyme table, BiG-SCAPE clustering statistics and GCF visualization
- **Novel BGCs**: BGC regions without KnownClusterBlast matches
- **KCB Hits**: Known cluster matches grouped by MIBiG entry

### Rarefaction Curve
- Shows GCF discovery saturation across sampled genomes
- Generated from BiG-SCAPE SQLite database (`{taxon}.db`)
- Displays total GCFs and saturation percentage
- Helps estimate diversity coverage and whether more sampling is needed

### GCF Visualization
- Shows representative BGCs for each Gene Cluster Family
- Includes gene arrows with functional annotations
- Color-coded by gene function (core biosynthetic, transport, regulatory, etc.)
- Links to antiSMASH results for detailed analysis (paths relative to `main_analysis_results/{taxon}/`)
- Displays KCB hit or "Potentially Novel" designation for each GCF representative
- Novel BGCs tab shows GCF family assignment when clustering is enabled

### GCF Biosynthetic Tree
- `GCF_BIOSYNTHETIC_TREE` runs **before** `VISUALIZE_RESULTS` — its PNG output is passed as `gcf_tree_png` input to create an explicit Nextflow data dependency
- The tree PNG is embedded as base64 in `bgc_report.html`, making the report self-contained
- Published copies also exist in `gcf_heatmap/` for standalone use
- `conf/conda.config` uses `withName: 'GCF_BIOSYNTHETIC_TREE'` for the conda environment

### BGC Distribution Analysis
- GCF × Genus heatmap showing BGC distribution across taxonomic groups
- Genus-specific GCFs table (potential taxon markers)
- Widespread GCFs table (found in 5+ genera, conserved or HGT)
- Uses GTDB-Tk taxonomy when available, falls back to NCBI taxonomy
- Phylogenetic tree files available in `results/gtdbtk_results/` for external viewers (iTOL, FigTree)

## Post-Pipeline BGC Analysis Scripts

These standalone scripts (in `scripts/`) perform additional analyses after the main pipeline completes. They operate on the BiG-SCAPE SQLite database and antiSMASH outputs.

### `scripts/bgc_pfam_tree.py` — Jaccard-distance NJ tree of BGCs

Builds a Neighbor-Joining tree based on Pfam domain presence/absence (Jaccard distance).

```bash
python scripts/bgc_pfam_tree.py \
    --db results/bigscape_results/Pantoea/Pantoea.db \
    --bgc_type phosphonate \
    --outdir results/bgc_trees/Pantoea
    [--family_id 2]    # Optional: restrict to a single GCF
```

Outputs: `_pfam_tree.nwk`, `_jaccard_distances.tsv`, `_domain_matrix.tsv`, `_metadata.json`

### `scripts/bgc_synteny_tree.py` — LCS-based gene-order tree

Builds a tree based on ordered domain sequences (one domain per CDS, sorted by genomic position). Uses normalized LCS distance. **Note:** Can be confused by strand orientation.

```bash
python scripts/bgc_synteny_tree.py \
    --db results/bigscape_results/Pantoea/Pantoea.db \
    --bgc_type phosphonate \
    --outdir results/bgc_trees/Pantoea/GCF2_synteny \
    [--family_id 2]
```

Outputs: `_synteny_tree.nwk`, `_domain_sequences.tsv`, `_lcs_distances.tsv`, `_metadata.json`

### `scripts/bgc_architecture_tree.py` — Architecture deduplication tree

Groups BGCs by exact domain multiset (orientation-independent), then builds a generalized Jaccard NJ tree of the unique architectures. Best for within-GCF comparison.

```bash
python scripts/bgc_architecture_tree.py \
    --db results/bigscape_results/Pantoea/Pantoea.db \
    --bgc_type phosphonate \
    --family_id 2 \
    --outdir results/bgc_trees/Pantoea/GCF2_arch
```

Outputs: NJ tree + five iTOL annotation files (count bar, domain binary, genus colorstrip, arch label, genome list).

Architecture labels: `arch_001_n138` (arch rank, count). GCF2 phosphonate → 205 BGCs → 21 unique architectures; arch_001 (n=138) is the dominant core.

### `scripts/bgc_coupling_annotation.py` — iTOL coupling enzyme colorstrip

Classifies each phosphonate BGC by the coupling enzyme acting on phosphonopyruvate (the branching step immediately downstream of PEP mutase). Reads antiSMASH JSON files for rich SMCOG and rule-based-cluster annotations, then outputs an iTOL DATASET_COLORSTRIP file.

```bash
python scripts/bgc_coupling_annotation.py \
    --antismash_dir results/antismash_results/Pantoea \
    --metadata results/bgc_trees/Pantoea/phosphonate_metadata.json \
    --outfile results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \
    --bgc_type phosphonate
```

**Coupling enzyme classes detected (Pantoea, n=1212 BGCs):**

| Class | Marker | Pathway | GCF | Count |
|-------|--------|---------|-----|-------|
| FrbC | SMCOG1271 (HMGL-like) | → phosphonomethylmalate → phosphinothricin-type | GCF-2/3 | 920 |
| Fe-ADH | Fe-ADH rule | → phosphonolactate (reductase route) | GCF-4/6 | 112 |
| TPP+NTP | TPP_enzyme_C + NTP_transf_3 rules | → phosphonolipid (CDP-pathway) | GCF-5 | 84 |
| Ppd | SMCOG1055 (ThDP-decarboxylase) | → 2-phosphonoacetaldehyde → 2-AEP | GCF-1/8 | 72 |
| PalB* | SMCOG1013 | → phosphonoalanine? | GCF-7 | 20 |
| Unknown | — | — | — | 4 |

**Key insights:**
- Classification maps almost perfectly onto BiG-SCAPE GCF families — coupling enzyme type is the primary determinant of GCF membership.
- The `Fe-ADH` rule-based marker (iron-containing alcohol dehydrogenase / 2-Hacid_dh_C) is antiSMASH's marker for the phosphonopyruvate reductase (→ phosphonolactate) pathway.
- GCF-5 (TPP+NTP) confirmed as **phosphonolipid BGCs**: Ppd-type ThDP enzyme + two NTP_transf_3 cytidylyltransferases + CDP-alcohol phosphatidyltransferases + Asn_synthase (CDP-phosphonate pathway). Well-annotated NCBI genomes explicitly label the ThDP enzyme as "phosphonopyruvate decarboxylase".
- AEP-pathway BGCs (GCF-1/8) use Ppd as coupling enzyme regardless of tailoring enzymes downstream.

**⚠️ PalB classification is pending correction:** The current script uses SMCOG1013 (Aminotran_3, fold type IV PLP) to detect PalB. However, PalB is an **AAT superfamily enzyme (fold type I PLP)** annotated as Aminotran_1_2 / PF00155 / SMCOG1019 — a completely different aminotransferase class. Additionally, coupling enzymes are not always adjacent to pepM in the BGC (the phosphonoalamide BGC architecture shows PalB far from pepM). The correct approach is protein sequence phylogenetic placement against characterized references (see "Coupling Enzyme Reference Trees" below).

**Note on PalA:** PalA (phosphonopyruvate hydrolase, a phosphonate degradation/resistance gene) does not confound the classification — all GCF types show clear biosynthetic markers.

### Coupling Enzyme Reference Trees

Protein sequence phylogenetic trees for each coupling enzyme class, with characterized MIBiG/literature reference sequences as anchors to place Pantoea BGC sequences.

**Reference sequences:** `results/bgc_trees/Pantoea/coupling_enzyme_trees/reference_coupling_enzymes.faa`

Seven characterized coupling enzymes covering all four reaction types:

| FASTA ID | Protein | Function | Source |
|----------|---------|----------|--------|
| `BGC0000904\|ABB90392\|FrbC` | FrbC | phosphonomethylmalate synthase | *Streptomyces rubellomurinus* (FR-900098) |
| `BGC0000897\|ACZ13457\|DhpF` | DhpF | phosphonopyruvate decarboxylase | *Streptomyces luridus* (Dehydrophos) |
| `BGC0000938\|ACG70832\|Fom2` | Fom2 | phosphonopyruvate decarboxylase | *Streptomyces fradiae* (Fosfomycin) |
| `BGC0000806\|AHL24480\|Ppd` | Ppd | phosphonopyruvate decarboxylase | *Glycomyces* sp. NRRL B-16210 |
| `Phosphonoalamide_BGC\|WP_051781701\|PnaA` | PnaA | phosphonopyruvate transaminase | *Streptomyces* sp. NRRL B-2790 |
| `Valinophos_BGC\|WP_063765859\|VlpB` | VlpB | phosphonopyruvate reductase | *Streptomyces durhamensis* NRRL B-3309 |
| `Pantaphos_BGC\|WP_013027159\|HvrC` | HvrC | phosphonomethylmalate synthase | *Pantoea ananatis* LMG 5342 |

**Plan (pending):** Extract coupling enzyme CDS sequences from Pantoea antiSMASH JSONs for each class, align with references using MAFFT/MUSCLE, build ML or NJ trees (FastTree/IQ-TREE), and use phylogenetic placement to:
1. Correctly classify GCF-7 (currently mislabeled "PalB" via SMCOG1013)
2. Confirm or reclassify the 4 "Unknown" BGCs
3. Identify true PalB-type (AAT superfamily, fold type I PLP) transaminase coupling enzymes in Pantoea

### `scripts/bgc_gcf_heatmap.py` — GCF × Species presence/absence heatmap

Generates a heatmap of GCF membership across organism groups, with a GTDB-Tk phylogenetic tree as column ordering and a Jaccard/complete-linkage row dendrogram matching BiG-SCAPE's clustering algorithm.

```bash
python scripts/bgc_gcf_heatmap.py \
    --metadata            results/bgc_trees/Pantoea/phosphonate_metadata.json \
    --coupling_annotation results/bgc_trees/Pantoea/phosphonate_itol_coupling.txt \
    --gtdbtk_tree         results/gtdbtk_results/Pantoea/gtdbtk_output/classify/gtdbtk.bac120.classify.tree.1.tree \
    --gtdbtk_summary      results/gtdbtk_results/Pantoea/gtdbtk_output/gtdbtk.bac120.summary.tsv \
    --outdir              results/bgc_trees/Pantoea
```

- **Data source**: Only region-level BGC records with GCF assignments at `cutoff=0.3` (303 phosphonate BGCs; sub-records like cand_cluster/protocluster are excluded)
- **True singletons**: Single-member GCFs (size=1), not unassigned records
- **Row dendrogram**: `scipy.spatial.distance.pdist(metric='jaccard')` + `linkage(method='complete')` — matches BiG-SCAPE's clustering algorithm
- **Column tree**: GTDB-Tk phylogenetic tree pruned to representative genomes per organism group, rendered as a cladogram
- **Outputs**: `gcf_species_heatmap.png` and `.svg`

### `scripts/bgc_itol_annotations.py` — iTOL annotation files

Generates iTOL annotation files from bgc_pfam_tree.py or bgc_synteny_tree.py outputs.

```bash
python scripts/bgc_itol_annotations.py \
    --treedir results/bgc_trees/Pantoea \
    --bgc_type phosphonate
```

Outputs: `_itol_gcf.txt` (color strip), `_itol_domains.txt` (binary), `_itol_domaincount.txt` (bar chart).

### Key Pfam accessions for phosphonate BGCs

Verified from antiSMASH clusterhmmer output on Pantoea phosphonate clusters:

| Pfam | Name | Function |
|------|------|----------|
| PF13714 | PEP_mutase | PEP mutase (pepM/aepX) — hallmark gene |
| PF00296 | HMGL-like (HEPD) | 2-hydroxyethylphosphonate dioxygenase |
| PF00682 | FrbC-like (PmmS) | Phosphonomethylmalate synthase (HMGL superfamily) |
| PF02775 | ThDP_C | Phosphonopyruvate decarboxylase |
| PF00266 | Aminotrans_V | 2-AEP transaminase |
| PF13649 | Radical_SAM | Radical C–P chemistry |

**Note on HMGL annotation:** AntiSMASH/BiG-SCAPE annotates phosphonomethylmalate synthase as `PF00682 (HMGL-like)` because it structurally belongs to the HMGL superfamily. The antiSMASH JSON provides richer context via `gene_functions: biosynthetic-additional (smcogs) SMCOG1271: 2-isopropylmalate synthase` and `sec_met_domain: HMGL-like`. BiG-SCAPE only stores the Pfam accession and bit score — no SMCOG or functional description.

### Data Sources

- BiG-SCAPE DB `hsp` table: Pfam accession + bit_score per CDS (populated by antiSMASH clusterhmmer)
- AntiSMASH JSON: richer annotations including `gene_functions`, `sec_met_domain` (SMCOG hits, TIGRFAM), and `product`
- Domain sequences in TSV come from the BiG-SCAPE DB (best Pfam hit per CDS, ordered by `nt_start`)

## Development Notes

### Configuration

- **Conda environments**: Defined centrally in `conf/conda.config` (not in individual modules)
- **Resource labels**: Defined in `conf/labels.config`, applied via process labels in modules
- **SLURM overrides**: Profile-specific adjustments in `nextflow.config`

### Utilities

- `Utils.sanitizeTaxon(name)`: Sanitize taxon for filesystem paths (removes special chars)
- `Utils.antismashParamsHash(params)`: Generate MD5 hash of antiSMASH parameters for reuse tracking
- `Utils.buildReusePath(params, projectDir, tool, taxon, subPath)`: Build absolute path for result reuse
- `Utils.isValidInput(input)`: Check if input is valid (not a placeholder)

### Module Guidelines

- Use `publishDir` for outputs, `storeDir` for database downloads
- Use appropriate labels: `process_low`, `process_medium`, `process_high`, `process_high_memory`
- Use `tolerant` label for per-genome processes where individual failures are acceptable
- antiSMASH uses `cache 'lenient'` for directory inputs
- COLLECT_VERSIONS searches `work/conda/` for installed tool versions
- BiG-SCAPE database (`bigscape_db`) is passed explicitly through pipeline for rarefaction curve generation

### Data Key Conventions

BGC regions are uniquely identified using `region_name` (e.g., "40.1" = record_index 40, region 1):

- **KCB lookup**: `(genome, region_name)` → KnownClusterBlast hit info
- **BGC-to-GCF mapping**: `(genome, region_name)` → GCF family assignment
- **JSON serialization**: `"genome|region_name"` format (e.g., `"Streptomyces_coelicolor_A32|40.1"`)

This avoids key collisions since `region` numbers are only unique within a record/contig, not within a genome. The `region_name` matches antiSMASH's naming convention directly.

Key files:
- `scripts/clustering/extract_gcf_representatives.py`: `load_kcb_lookup()`, `build_record_index_map()`, `extract_genome_gcf_mapping()`
- `scripts/analysis/tabulate_regions.py`: Creates `region_name` column in tabulation

### antiSMASH Link Paths

Links from `bgc_report.html` to antiSMASH results use paths relative to `main_analysis_results/{taxon}/`:
```
../../antismash_results/{taxon}/{genome}/index.html#r{record_index}c{region_number}
```
This is set in `scripts/clustering/extract_gcf_representatives.py` (`antismash_link`). If the report location changes, this depth must be updated accordingly.

### Known Issues

- **Duplicate gene names**: Some NCBI genomes have duplicate CDS feature names (e.g., `sapC`), causing antiSMASH to fail with "multiple CDS features have the same name"
- **DIAMOND memory errors**: `malloc(): corrupted top size` errors during ClusterBlast indicate memory issues; try increasing memory allocation or reducing concurrent jobs
- **NCBI dehydrated download corruption**: The `--dehydrated` download mode can produce null-filled files due to network timeouts. The module includes validation with `sync` + retry logic, but if corruption persists, delete the cached work directory and re-run
- **pyhmmer/BiG-SCAPE compatibility**: pyhmmer 0.12+ changed `profile.accession` from bytes to str, breaking BiG-SCAPE. The module pins `pyhmmer<0.11`
- **GTDB-Tk duplicate taxon labels**: GTDB-Tk normalizes genome names case-insensitively. If two genomes have names differing only in case (e.g., `MDCuke` vs `MDcuke`), GTDB-Tk will fail with `NewickReaderDuplicateTaxonError`. The `create_name_map.py` script now handles this by tracking names case-insensitively and adding numeric suffixes to duplicates

## Troubleshooting

```bash
nextflow clean -f -k              # Clear cache if -resume fails
rm -rf work/                      # Remove intermediate files (keeps databases)
du -sh work/                      # Check work dir size
```

- **GTDB-Tk OOM**: Requires 56-64GB RAM; keep `--pplacer_cpus 1`
- **Conda env conflicts**: Each tool has its own environment; don't mix in COLLECT_VERSIONS
- **NCBI download corruption**: If GENBANK_TO_FASTA fails with "No sequences found", check for null-filled files:
  ```bash
  # Find corrupted files (first bytes are null)
  for f in work/*/ncbi_dataset/data/*/genomic.gbff; do
    [ -z "$(head -c 10 "$f" | tr -d '\0')" ] && echo "Corrupted: $f"
  done
  # Fix: delete cached NCBI download and re-run
  grep "NCBI_DATASETS_DOWNLOAD" .nextflow.log | grep "workDir" | tail -1  # Find work dir
  rm -rf work/XX/XXXXXX  # Delete the cached directory
  nextflow run main.nf -resume
  ```
