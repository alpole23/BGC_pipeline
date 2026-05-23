# ClusterQuest Parameters Reference

This page documents all configurable and non-configurable parameters in ClusterQuest. Parameters are set via the command line (`--param value`) or by editing `nextflow.config` directly.

---

## Configurable Parameters

### Global

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--taxon` | `"Erwiniaceae"` | string | **Required.** NCBI taxon name (species, genus, family, order, etc.) |
| `--workflow` | `"full"` | string | Pipeline mode: `download`, `bgc_analysis`, or `full` |
| `--outdir` | `"results"` | string | Output directory for all results |

**Workflow modes:**
- `download` — Download and prepare genomes from NCBI only
- `bgc_analysis` — Run BGC analysis on pre-downloaded genomes
- `full` — End-to-end: download + BGC analysis (recommended)

---

### Genome Download / Input

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--input_genomes` | `null` | path | Path to pre-downloaded genomes directory. Only used with `--workflow bgc_analysis` to skip NCBI download |

---

### antiSMASH (BGC Detection)

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--antismash_minimal` | `false` | bool | Minimal mode: skips domain analysis for faster runs. Disables most `--antismash_*` options below |
| `--antismash_cb_general` | `false` | bool | ClusterBlast: compare detected BGCs against the full antiSMASH database |
| `--antismash_cc_mibig` | `false` | bool | ClusterCompare: advanced scoring against MIBiG (more sensitive than KnownClusterBlast) |
| `--antismash_smcog_trees` | `true` | bool | Generate phylogenetic trees for BGC core biosynthetic genes |
| `--reuse_antismash_from` | `null` | string | Taxon name to reuse antiSMASH results from (see [Cross-Taxon Reuse](#cross-taxon-result-reuse)) |

> **Note:** Detection is hardcoded to phosphonate BGCs only (`--hmmdetection-limit-to-rule-names phosphonate`). KnownClusterBlast (`--cb-knownclusters`), `clusterhmmer`, and `tigrfam` are always enabled and cannot be turned off.

---

### Region Analysis

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--run_analysis` | `true` | bool | Enable BGC counting, tabulation, and HTML report generation |
| `--count_per_contig` | `false` | bool | Count BGC regions per contig (`true`) or aggregate per genome (`false`) |
| `--split_hybrids` | `false` | bool | Split hybrid BGC types into components (e.g., `T1PKS-NRPS` → `T1PKS` + `NRPS`) |

---

### Clustering

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--clustering` | `"bigscape"` | string | Clustering tool: `none` or `bigscape` |

#### BiG-SCAPE Options
Active when `--clustering bigscape`.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--bigscape_cutoffs` | `"0.30"` | string | GCF distance threshold(s); comma-separated for multiple (e.g., `"0.30,0.40"`) |
| `--bigscape_alignment_mode` | `"auto"` | string | Alignment mode: `auto`, `global`, or `glocal` |
| `--bigscape_mibig_version` | `""` | string | Include MIBiG reference clusters in network (e.g., `"3.1"`). Empty string excludes MIBiG |
| `--bigscape_classify` | `"category"` | string | BGC classification scheme: `""` (none), `"category"`, `"class"`, or `"legacy"` |
| `--bigscape_include_singletons` | `true` | bool | Include unclustered BGCs (singletons) in GCF output |
| `--bigscape_mix` | `false` | bool | Allow BGCs of different classes to cluster together in the same GCF |

---

### Phylogeny (GTDB-Tk)

> **Requirements:** ~140 GB disk space, 56–64 GB RAM. Keep `--gtdbtk_pplacer_cpus 1` to avoid out-of-memory errors.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `--run_gtdbtk` | `true` | bool | Enable phylogenetic placement with GTDB-Tk |
| `--gtdbtk_bgc_genomes_only` | `true` | bool | Only analyze genomes that have at least one detected BGC (reduces runtime) |
| `--gtdbtk_cpus` | `8` | int | CPUs for GTDB-Tk identify and align steps |
| `--gtdbtk_pplacer_cpus` | `1` | int | CPUs for pplacer tree placement (memory-bound; keep at 1 unless high-RAM system) |
| `--gtdbtk_min_perc_aa` | `10` | int | Minimum percentage of amino acids in the multiple-sequence alignment to include a genome |
| `--gtdbtk_outgroup` | `null` | string | Outgroup pattern for tree rooting (e.g., `"g__Escherichia"`). `null` = no rooting |
| `--reuse_gtdbtk_from` | `null` | string | Taxon name to reuse GTDB-Tk results from (see [Cross-Taxon Reuse](#cross-taxon-result-reuse)) |

---

## Cross-Taxon Result Reuse

ClusterQuest can reuse antiSMASH and GTDB-Tk results from a previous run to avoid redundant computation when analyzing related taxa.

```bash
# First: run on a broad taxon
nextflow run main.nf --taxon "Erwiniaceae"

# Later: reuse results for a subset taxon
nextflow run main.nf --taxon "Pantoea" \
  --reuse_antismash_from "Erwiniaceae" \
  --reuse_gtdbtk_from "Erwiniaceae"
```

**Reuse behavior:**

| Feature | Direction | Behavior |
|---------|-----------|----------|
| antiSMASH reuse | Both directions (superset→subset or subset→superset) | Per-genome: each genome is checked individually. Matching genomes are copied; missing ones are run fresh. Partial reuse works. |
| GTDB-Tk reuse | Superset → subset only | Requires ALL current genomes to exist in the source results. If any genome is missing, falls back to a full GTDB-Tk run. |

---

## Non-Configurable Settings

These settings are hard-coded in the pipeline and cannot be changed via command-line parameters. They can only be modified by editing the relevant configuration file directly.

### Fixed antiSMASH Behavior

| Behavior | Source | Reason |
|----------|--------|--------|
| Phosphonate detection only | `modules/analysis/antismash.nf` | Hardcoded to `--hmmdetection-limit-to-rule-names phosphonate` |
| `--cb-knownclusters` always enabled | `modules/analysis/antismash.nf` | Always compare against MIBiG |
| `clusterhmmer` always enabled | `modules/analysis/antismash.nf` | Ensures consistent domain annotation across all analyses |
| `tigrfam` always enabled | `modules/analysis/antismash.nf` | Required for reliable BGC gene family classification |

### Resource Labels (`conf/labels.config`)

Process resource allocations are controlled by labels. These apply to all processes unless overridden by a profile (e.g., `slurm`).

| Label | CPUs | Memory | Time | Used by |
|-------|------|--------|------|---------|
| `process_local` | 1 | 1 GB | — | Download/setup steps (runs on head node) |
| `process_low` | 1 | 2 GB | 1h | Python scripts, tabulation, visualization |
| `process_medium` | 4 | 8 GB | 2h | antiSMASH (per-genome BGC detection) |
| `process_high` | 8 | 32 GB | 8h | BiG-SCAPE clustering |
| `process_high_memory` | 8 | 48 GB* | 24h | GTDB-Tk pplacer |

*The SLURM profile overrides `process_high_memory` to 128 GB.

### SLURM Profile Overrides (`nextflow.config`)

When running with `-profile slurm`, the following process-specific overrides apply on top of the resource labels:

| Process | Override |
|---------|----------|
| `BIGSCAPE` | Memory: 64 GB |
| `VISUALIZE_RESULTS` | Memory: 16 GB, Time: 4h |
| `DOWNLOAD_GTDBTK_DB` | Time: 12h |
| `NCBI_DATASETS_DOWNLOAD` | Time: 4h |
| `process_high_memory` label | Queue: `highmem`, Memory: 128 GB |

The SLURM executor is configured with:
- Queue size: 200 concurrent jobs
- Submit rate limit: 20 jobs/minute

### Pipeline Reporting (`nextflow.config`)

Nextflow automatically generates execution reports in `results/pipeline_info/`:

| Output | File | Content |
|--------|------|---------|
| Trace | `pipeline_trace.tsv` | Per-task CPU, memory, runtime, I/O stats |
| Report | `pipeline_report.html` | Aggregated resource usage report |
| Timeline | `pipeline_timeline.html` | Gantt chart of process execution |

These are always enabled and cannot be disabled via parameters (edit `nextflow.config` to change).

### Conda Configuration

Conda environments are always enabled and use channels in this priority order:
1. `conda-forge`
2. `bioconda`
3. `defaults`

Environment definitions are in `conf/conda.config`.

---

## Example Commands

```bash
# Minimal run (just BGC detection, no clustering or phylogeny)
nextflow run main.nf --taxon "Pantoea" --clustering none --run_gtdbtk false

# Full run with clustering and MIBiG
nextflow run main.nf --taxon "Streptomyces" \
  --bigscape_mibig_version "3.1" \
  --bigscape_cutoffs "0.30,0.40"

# Resume a previous run
nextflow run main.nf -resume

# HPC execution
nextflow run main.nf -profile slurm --taxon "Erwiniaceae"
```
