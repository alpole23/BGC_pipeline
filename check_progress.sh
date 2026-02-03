#!/bin/bash
# =============================================================================
# ClusterQuest Pipeline Progress Monitor
# =============================================================================
#
# Usage:
#   ./check_progress.sh              # Check progress with defaults
#   ./check_progress.sh -w           # Watch mode (auto-refresh every 10s)
#   ./check_progress.sh -d results2  # Use different output directory
#
# =============================================================================

set -o pipefail

# ----------------Configuration-------------------
OUTDIR="results"
WATCH_MODE=false
WATCH_INTERVAL=10

# Parse arguments
while getopts "wd:h" opt; do
    case $opt in
        w) WATCH_MODE=true ;;
        d) OUTDIR="$OPTARG" ;;
        h) echo "Usage: $0 [-w] [-d outdir]"; exit 0 ;;
        *) echo "Usage: $0 [-w] [-d outdir]"; exit 1 ;;
    esac
done

TRACE_FILE="${OUTDIR}/pipeline_info/pipeline_trace.tsv"
LOG_FILE=".nextflow.log"

# ----------------Functions-----------------------
print_header() {
    echo "╔═════════════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                           ClusterQuest Pipeline Progress                                 ║"
    echo "╠═════════════════════════════════════════════════════════════════════════════════════════╣"
    printf "║  %-86s ║\n" "$(date)"
    echo "╚═════════════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
}

print_summary() {
    echo "┌─────────────────────────────────────────────────────────────────────────────────────────┐"
    echo "│ Summary                                                                                 │"
    echo "├─────────────────────────────────────────────────────────────────────────────────────────┤"

    if [[ -f "$LOG_FILE" ]]; then
        SUBMITTED=$(grep -c 'Submitted process' "$LOG_FILE" 2>/dev/null || echo 0)
        COMPLETED=$(grep -c 'COMPLETED' "$LOG_FILE" 2>/dev/null || echo 0)
        FAILED=$(grep -c 'FAILED' "$LOG_FILE" 2>/dev/null || echo 0)
        CACHED=$(grep -c 'CACHED' "$LOG_FILE" 2>/dev/null || echo 0)

        printf "│  Submitted: %-10s  Completed: %-10s  Cached: %-10s  Failed: %-10s │\n" \
            "$SUBMITTED" "$COMPLETED" "$CACHED" "$FAILED"
    else
        printf "│  %-86s │\n" "No log file found (.nextflow.log)"
    fi
    echo "└─────────────────────────────────────────────────────────────────────────────────────────┘"
    echo ""
}

print_progress_bars() {
    echo "┌─────────────────────────────────────────────────────────────────────────────────────────┐"
    echo "│ Progress by Process                                                                     │"
    echo "├─────────────────────────────────────────────────────────────────────────────────────────┤"

    if [[ ! -f "$LOG_FILE" ]]; then
        printf "│  %-86s │\n" "Waiting for pipeline to start..."
        echo "└─────────────────────────────────────────────────────────────────────────────────────────┘"
        return
    fi

    # Process log to get counts
    {
        # Get submitted counts
        grep 'Submitted process' "$LOG_FILE" 2>/dev/null | \
            sed 's/.*> //;s/ (.*//' | sort | uniq -c | \
            awk '{print "S|" $2 "|" $1}'

        # Get completed counts from trace
        if [[ -f "$TRACE_FILE" ]]; then
            tail -n +2 "$TRACE_FILE" 2>/dev/null | \
                awk -F'\t' '$4=="COMPLETED" {print $3}' | \
                sed 's/ (.*//' | sort | uniq -c | \
                awk '{print "C|" $2 "|" $1}'
        fi
    } | awk -F'|' '
    {
        proc = $2
        count = $3
        if ($1 == "S") submitted[proc] = count
        if ($1 == "C") completed[proc] = count
    }
    END {
        for (proc in submitted) {
            s = submitted[proc]
            c = (proc in completed) ? completed[proc] : 0
            total = (c > s) ? c : s
            pct = (total > 0) ? int(c * 100 / total) : 0
            if (pct > 100) pct = 100

            # Build progress bar
            filled = int(pct * 25 / 100)
            bar = ""
            for (i = 0; i < filled; i++) bar = bar "█"
            for (i = filled; i < 25; i++) bar = bar "░"

            # Show full process name (up to 42 chars to fit box with counts)
            display_proc = proc
            if (length(proc) > 42) {
                display_proc = substr(proc, 1, 39) "..."
            }

            # Format: │  name [bar] pct% (n/t)  │  = 91 chars total
            printf "│  %-42s [%s] %3d%% (%d/%d)\n", display_proc, bar, pct, c, total
        }
    }' | sort

    echo "└─────────────────────────────────────────────────────────────────────────────────────────┘"
    echo ""
}

print_running() {
    echo "┌─────────────────────────────────────────────────────────────────────────────────────────┐"
    echo "│ Currently Running                                                                       │"
    echo "├─────────────────────────────────────────────────────────────────────────────────────────┤"

    if [[ -f "$TRACE_FILE" ]]; then
        RUNNING=$(tail -n +2 "$TRACE_FILE" 2>/dev/null | awk -F'\t' '$4=="-" || $4=="RUNNING" {print $3}' | head -5)
        if [[ -n "$RUNNING" ]]; then
            echo "$RUNNING" | while read -r task; do
                printf "│  %-86s │\n" "$task"
            done
        else
            printf "│  %-86s │\n" "None (pipeline may be between tasks or waiting for resources)"
        fi
    else
        printf "│  %-86s │\n" "Trace file not available yet"
    fi

    echo "└─────────────────────────────────────────────────────────────────────────────────────────┘"
    echo ""
}

print_recent() {
    echo "┌─────────────────────────────────────────────────────────────────────────────────────────┐"
    echo "│ Recent Activity                                                                         │"
    echo "├─────────────────────────────────────────────────────────────────────────────────────────┤"

    if [[ -f "$LOG_FILE" ]]; then
        tail -50 "$LOG_FILE" 2>/dev/null | grep -E "(Submitted|COMPLETED|FAILED|CACHED)" | tail -5 | \
        while read -r line; do
            # Extract just the process name and status
            short=$(echo "$line" | sed 's/.*\[\([^]]*\)\].*/\1/' | cut -c1-83)
            printf "│  %-86s │\n" "$short"
        done
    else
        printf "│  %-86s │\n" "No activity yet"
    fi

    echo "└─────────────────────────────────────────────────────────────────────────────────────────┘"
    echo ""
}

print_errors() {
    if [[ -f "$LOG_FILE" ]] && grep -q "FAILED\|ERROR" "$LOG_FILE" 2>/dev/null; then
        echo "┌─────────────────────────────────────────────────────────────────────────────────────────┐"
        echo "│ Errors                                                                                  │"
        echo "├─────────────────────────────────────────────────────────────────────────────────────────┤"

        grep -E "FAILED|ERROR" "$LOG_FILE" 2>/dev/null | tail -3 | \
        while read -r line; do
            short=$(echo "$line" | cut -c1-83)
            printf "│  %-86s │\n" "$short"
        done

        echo "└─────────────────────────────────────────────────────────────────────────────────────────┘"
        echo ""
    fi
}

# ----------------Main----------------------------
show_progress() {
    clear 2>/dev/null || true
    print_header
    print_summary
    print_progress_bars
    print_running
    print_recent
    print_errors
}

if $WATCH_MODE; then
    echo "Watch mode enabled. Press Ctrl+C to exit."
    while true; do
        show_progress
        echo "Refreshing in ${WATCH_INTERVAL}s... (Ctrl+C to exit)"
        sleep $WATCH_INTERVAL
    done
else
    show_progress
fi
