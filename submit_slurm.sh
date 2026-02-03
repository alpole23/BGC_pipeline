#!/bin/bash
# =============================================================================
# ClusterQuest SLURM Submission Script
# =============================================================================
#
# Usage:
#   sbatch submit_slurm.sh                    # Run with defaults
#   sbatch submit_slurm.sh "Streptomyces"     # Specify taxon
#   sbatch submit_slurm.sh "Pantoea" bigscape # Specify taxon and clustering
#
# =============================================================================

# ----------------SLURM Parameters----------------
#SBATCH -p normal                    # Partition (adjust for your cluster)
#SBATCH -n 2                         # CPUs for Nextflow head process
#SBATCH --mem=8g                     # Memory for head job
#SBATCH -N 1                         # Single node for head job
#SBATCH --time=48:00:00              # Max wall time
#SBATCH --mail-type=END,FAIL         # Email on completion or failure
#SBATCH -J clusterquest               # Job name
#SBATCH -o clusterquest_%j.out        # Standard output
#SBATCH -e clusterquest_%j.err        # Standard error

# ----------------Configuration-------------------
# Modify these defaults or pass as arguments

TAXON="${1:-Pantoea}"                # First argument or default
CLUSTERING="${2:-bigscape}"          # Second argument or default (none/bigscape/bigslice/both)
OUTDIR="results"
RESUME="-resume"                     # Set to "" to start fresh

# ----------------Environment Setup---------------
# Uncomment and modify for your cluster

# Option 1: Load modules (cluster-specific)
# module load nextflow/24.10.4
# module load anaconda3

# Option 2: Use conda environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate nextflow

# Nextflow configuration
export NXF_HOME="${PWD}/.nextflow"
export NXF_DISABLE_CHECK_LATEST=true
export NXF_ANSI_LOG=false            # Disable ANSI colors in log files

# ----------------Validation----------------------
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: nextflow not found. Load the module or activate conda environment."
    exit 1
fi

# ----------------Run Pipeline--------------------
echo "=============================================="
echo "ClusterQuest Pipeline"
echo "=============================================="
echo "Start time:  $(date)"
echo "Taxon:       ${TAXON}"
echo "Clustering:  ${CLUSTERING}"
echo "Output:      ${OUTDIR}"
echo "Job ID:      ${SLURM_JOB_ID:-local}"
echo "Node:        ${SLURMD_NODENAME:-$(hostname)}"
echo "=============================================="
echo ""

nextflow run main.nf \
    -profile slurm \
    --taxon "${TAXON}" \
    --clustering "${CLUSTERING}" \
    --outdir "${OUTDIR}" \
    ${RESUME}

EXIT_CODE=$?

# ----------------Summary-------------------------
echo ""
echo "=============================================="
echo "Pipeline Complete"
echo "=============================================="
echo "End time:    $(date)"
echo "Exit code:   ${EXIT_CODE}"

if [ -f "${OUTDIR}/pipeline_info/pipeline_trace.tsv" ]; then
    TOTAL_TASKS=$(tail -n +2 "${OUTDIR}/pipeline_info/pipeline_trace.tsv" | wc -l)
    COMPLETED=$(tail -n +2 "${OUTDIR}/pipeline_info/pipeline_trace.tsv" | grep -c "COMPLETED")
    FAILED=$(tail -n +2 "${OUTDIR}/pipeline_info/pipeline_trace.tsv" | grep -c "FAILED")
    echo "Tasks:       ${COMPLETED}/${TOTAL_TASKS} completed, ${FAILED} failed"
fi

if [ -f "${OUTDIR}/main_analysis_results/$(echo ${TAXON} | tr ' ' '_')/main_data_visualization/bgc_report.html" ]; then
    echo "Report:      ${OUTDIR}/main_analysis_results/*/main_data_visualization/bgc_report.html"
fi

echo "=============================================="

exit ${EXIT_CODE}
