#!/usr/bin/env python3
"""
Resource Usage Visualization Script

Parses Nextflow trace files and generates an interactive HTML report
showing memory, CPU, and runtime statistics per process.

Usage:
    python visualize_resources.py --trace pipeline_trace.tsv --outdir results/
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import timedelta
import re


def parse_duration(duration_str):
    """Parse Nextflow duration string to seconds."""
    if not duration_str or duration_str == '-':
        return 0

    # Handle format like "1h 23m 45s" or "23m 45s" or "45s" or "1.5s"
    total_seconds = 0

    # Try parsing as milliseconds first (format: 123ms)
    if 'ms' in duration_str:
        match = re.search(r'([\d.]+)ms', duration_str)
        if match:
            return float(match.group(1)) / 1000

    # Parse hours
    match = re.search(r'(\d+)h', duration_str)
    if match:
        total_seconds += int(match.group(1)) * 3600

    # Parse minutes
    match = re.search(r'(\d+)m(?!s)', duration_str)
    if match:
        total_seconds += int(match.group(1)) * 60

    # Parse seconds
    match = re.search(r'([\d.]+)s', duration_str)
    if match:
        total_seconds += float(match.group(1))

    return total_seconds if total_seconds > 0 else 0


def parse_memory(mem_str):
    """Parse memory string to bytes."""
    if not mem_str or mem_str == '-':
        return 0

    mem_str = str(mem_str).strip().upper()

    # Handle numeric values (already in bytes)
    try:
        return int(float(mem_str))
    except ValueError:
        pass

    # Parse with units
    units = {
        'B': 1,
        'KB': 1024,
        'MB': 1024**2,
        'GB': 1024**3,
        'TB': 1024**4,
        'K': 1024,
        'M': 1024**2,
        'G': 1024**3,
        'T': 1024**4,
    }

    match = re.match(r'([\d.]+)\s*([KMGT]?B?)', mem_str)
    if match:
        value = float(match.group(1))
        unit = match.group(2) or 'B'
        return int(value * units.get(unit, 1))

    return 0


def format_bytes(bytes_val):
    """Format bytes to human-readable string."""
    if bytes_val == 0:
        return "0 B"

    units = ['B', 'KB', 'MB', 'GB', 'TB']
    unit_idx = 0
    val = float(bytes_val)

    while val >= 1024 and unit_idx < len(units) - 1:
        val /= 1024
        unit_idx += 1

    return f"{val:.2f} {units[unit_idx]}"


def format_duration(seconds):
    """Format seconds to human-readable duration."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        mins = int(seconds // 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds // 3600)
        mins = int((seconds % 3600) // 60)
        return f"{hours}h {mins}m"


def parse_trace_file(trace_path):
    """Parse Nextflow trace TSV file."""
    tasks = []

    with open(trace_path, 'r') as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')

            if header is None:
                header = parts
                continue

            task = dict(zip(header, parts))
            tasks.append(task)

    return tasks


def aggregate_by_process(tasks):
    """Aggregate task statistics by process name."""
    processes = {}

    for task in tasks:
        name = task.get('name', 'unknown')
        # Extract base process name (remove parameters in parentheses)
        base_name = re.sub(r'\s*\([^)]*\)', '', name).strip()

        if base_name not in processes:
            processes[base_name] = {
                'name': base_name,
                'count': 0,
                'completed': 0,
                'failed': 0,
                'total_runtime': 0,
                'peak_memory': 0,
                'total_memory': 0,
                'peak_cpu': 0,
                'total_cpu': 0,
                'memory_requested': 0,
                'tasks': []
            }

        proc = processes[base_name]
        proc['count'] += 1

        status = task.get('status', '')
        if status == 'COMPLETED':
            proc['completed'] += 1
        elif status in ('FAILED', 'ABORTED'):
            proc['failed'] += 1

        # Parse metrics
        runtime = parse_duration(task.get('realtime', '0'))
        peak_rss = parse_memory(task.get('peak_rss', '0'))
        cpu_pct = float(task.get('%cpu', '0').replace('%', '') or 0)
        mem_requested = parse_memory(task.get('memory', '0'))

        proc['total_runtime'] += runtime
        proc['total_memory'] += peak_rss
        proc['total_cpu'] += cpu_pct
        proc['peak_memory'] = max(proc['peak_memory'], peak_rss)
        proc['peak_cpu'] = max(proc['peak_cpu'], cpu_pct)
        proc['memory_requested'] = max(proc['memory_requested'], mem_requested)

        proc['tasks'].append({
            'name': name,
            'status': status,
            'runtime': runtime,
            'peak_rss': peak_rss,
            'cpu_pct': cpu_pct,
            'memory_requested': mem_requested
        })

    # Calculate averages
    for proc in processes.values():
        if proc['count'] > 0:
            proc['avg_runtime'] = proc['total_runtime'] / proc['count']
            proc['avg_memory'] = proc['total_memory'] / proc['count']
            proc['avg_cpu'] = proc['total_cpu'] / proc['count']
            if proc['memory_requested'] > 0:
                proc['memory_efficiency'] = (proc['peak_memory'] / proc['memory_requested']) * 100
            else:
                proc['memory_efficiency'] = 0

    return processes


def generate_html_report(processes, output_path, trace_path):
    """Generate interactive HTML report."""

    # Sort processes by total runtime
    sorted_procs = sorted(processes.values(), key=lambda x: x['total_runtime'], reverse=True)

    # Calculate totals
    total_tasks = sum(p['count'] for p in sorted_procs)
    total_completed = sum(p['completed'] for p in sorted_procs)
    total_failed = sum(p['failed'] for p in sorted_procs)
    total_runtime = sum(p['total_runtime'] for p in sorted_procs)
    peak_memory_overall = max((p['peak_memory'] for p in sorted_procs), default=0)

    # Prepare data for charts
    process_names = [p['name'] for p in sorted_procs]
    runtime_data = [p['total_runtime'] / 60 for p in sorted_procs]  # Convert to minutes
    memory_data = [p['peak_memory'] / (1024**3) for p in sorted_procs]  # Convert to GB
    task_counts = [p['count'] for p in sorted_procs]
    efficiency_data = [p.get('memory_efficiency', 0) for p in sorted_procs]

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pipeline Resource Usage Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1 {{
            text-align: center;
            color: #2c3e50;
            margin-bottom: 30px;
            padding: 20px;
            background: white;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .summary-card {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .summary-card h3 {{
            color: #7f8c8d;
            font-size: 0.9em;
            margin-bottom: 10px;
        }}
        .summary-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
        }}
        .summary-card .value.success {{
            color: #27ae60;
        }}
        .summary-card .value.danger {{
            color: #e74c3c;
        }}
        .tabs {{
            display: flex;
            gap: 10px;
            margin-bottom: 20px;
            flex-wrap: wrap;
        }}
        .tab-btn {{
            padding: 10px 20px;
            border: none;
            background: white;
            cursor: pointer;
            border-radius: 5px;
            font-size: 1em;
            transition: all 0.3s;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .tab-btn:hover {{
            background: #ecf0f1;
        }}
        .tab-btn.active {{
            background: #3498db;
            color: white;
        }}
        .tab-content {{
            display: none;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .tab-content.active {{
            display: block;
        }}
        .chart-container {{
            position: relative;
            height: 400px;
            margin-bottom: 30px;
        }}
        .chart-row {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 30px;
            margin-bottom: 30px;
        }}
        @media (max-width: 900px) {{
            .chart-row {{
                grid-template-columns: 1fr;
            }}
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }}
        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #2c3e50;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .status-completed {{
            color: #27ae60;
            font-weight: bold;
        }}
        .status-failed {{
            color: #e74c3c;
            font-weight: bold;
        }}
        .progress-bar {{
            width: 100%;
            height: 20px;
            background: #ecf0f1;
            border-radius: 10px;
            overflow: hidden;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #3498db, #2ecc71);
            border-radius: 10px;
            transition: width 0.3s;
        }}
        .efficiency-low {{
            background: linear-gradient(90deg, #e74c3c, #e67e22) !important;
        }}
        .efficiency-medium {{
            background: linear-gradient(90deg, #f1c40f, #2ecc71) !important;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Pipeline Resource Usage Report</h1>

        <div class="summary-grid">
            <div class="summary-card">
                <h3>Total Tasks</h3>
                <div class="value">{total_tasks}</div>
            </div>
            <div class="summary-card">
                <h3>Completed</h3>
                <div class="value success">{total_completed}</div>
            </div>
            <div class="summary-card">
                <h3>Failed</h3>
                <div class="value {"danger" if total_failed > 0 else ""}">{total_failed}</div>
            </div>
            <div class="summary-card">
                <h3>Total Runtime</h3>
                <div class="value">{format_duration(total_runtime)}</div>
            </div>
            <div class="summary-card">
                <h3>Peak Memory</h3>
                <div class="value">{format_bytes(peak_memory_overall)}</div>
            </div>
            <div class="summary-card">
                <h3>Processes</h3>
                <div class="value">{len(sorted_procs)}</div>
            </div>
        </div>

        <div class="tabs">
            <button class="tab-btn active" onclick="showTab('overview')">Overview</button>
            <button class="tab-btn" onclick="showTab('memory')">Memory Usage</button>
            <button class="tab-btn" onclick="showTab('runtime')">Runtime</button>
            <button class="tab-btn" onclick="showTab('efficiency')">Efficiency</button>
            <button class="tab-btn" onclick="showTab('details')">Process Details</button>
        </div>

        <div id="overview" class="tab-content active">
            <h2>Resource Overview</h2>
            <div class="chart-row">
                <div class="chart-container">
                    <canvas id="runtimeChart"></canvas>
                </div>
                <div class="chart-container">
                    <canvas id="memoryChart"></canvas>
                </div>
            </div>
            <div class="chart-container">
                <canvas id="taskCountChart"></canvas>
            </div>
        </div>

        <div id="memory" class="tab-content">
            <h2>Memory Usage by Process</h2>
            <div class="chart-container" style="height: 500px;">
                <canvas id="memoryDetailChart"></canvas>
            </div>
            <table>
                <thead>
                    <tr>
                        <th>Process</th>
                        <th>Peak Memory</th>
                        <th>Avg Memory</th>
                        <th>Requested</th>
                        <th>Efficiency</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(f'''
                    <tr>
                        <td>{p['name']}</td>
                        <td>{format_bytes(p['peak_memory'])}</td>
                        <td>{format_bytes(p['avg_memory'])}</td>
                        <td>{format_bytes(p['memory_requested'])}</td>
                        <td>
                            <div class="progress-bar">
                                <div class="progress-fill {"efficiency-low" if p.get('memory_efficiency', 0) < 30 else "efficiency-medium" if p.get('memory_efficiency', 0) < 70 else ""}"
                                     style="width: {min(p.get('memory_efficiency', 0), 100):.0f}%"></div>
                            </div>
                            {p.get('memory_efficiency', 0):.1f}%
                        </td>
                    </tr>
                    ''' for p in sorted_procs if p['count'] > 0)}
                </tbody>
            </table>
        </div>

        <div id="runtime" class="tab-content">
            <h2>Runtime by Process</h2>
            <div class="chart-container" style="height: 500px;">
                <canvas id="runtimeDetailChart"></canvas>
            </div>
            <table>
                <thead>
                    <tr>
                        <th>Process</th>
                        <th>Total Runtime</th>
                        <th>Avg Runtime</th>
                        <th>Tasks</th>
                        <th>Avg CPU %</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(f'''
                    <tr>
                        <td>{p['name']}</td>
                        <td>{format_duration(p['total_runtime'])}</td>
                        <td>{format_duration(p['avg_runtime'])}</td>
                        <td>{p['count']}</td>
                        <td>{p['avg_cpu']:.1f}%</td>
                    </tr>
                    ''' for p in sorted_procs if p['count'] > 0)}
                </tbody>
            </table>
        </div>

        <div id="efficiency" class="tab-content">
            <h2>Resource Efficiency</h2>
            <p style="margin-bottom: 20px; color: #7f8c8d;">
                Memory efficiency shows how much of the requested memory was actually used (peak usage / requested).
                Low efficiency suggests over-allocation; consider reducing memory requests.
            </p>
            <div class="chart-container" style="height: 500px;">
                <canvas id="efficiencyChart"></canvas>
            </div>
        </div>

        <div id="details" class="tab-content">
            <h2>Process Details</h2>
            <table>
                <thead>
                    <tr>
                        <th>Process</th>
                        <th>Tasks</th>
                        <th>Completed</th>
                        <th>Failed</th>
                        <th>Total Runtime</th>
                        <th>Peak Memory</th>
                        <th>Peak CPU</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(f'''
                    <tr>
                        <td>{p['name']}</td>
                        <td>{p['count']}</td>
                        <td class="status-completed">{p['completed']}</td>
                        <td class="{"status-failed" if p['failed'] > 0 else ""}">{p['failed']}</td>
                        <td>{format_duration(p['total_runtime'])}</td>
                        <td>{format_bytes(p['peak_memory'])}</td>
                        <td>{p['peak_cpu']:.1f}%</td>
                    </tr>
                    ''' for p in sorted_procs)}
                </tbody>
            </table>
        </div>

        <div class="footer">
            <p>Generated from: {trace_path}</p>
            <p>BGC Analysis Pipeline - Resource Usage Report</p>
        </div>
    </div>

    <script>
        const processNames = {json.dumps(process_names)};
        const runtimeData = {json.dumps(runtime_data)};
        const memoryData = {json.dumps(memory_data)};
        const taskCounts = {json.dumps(task_counts)};
        const efficiencyData = {json.dumps(efficiency_data)};

        const colors = [
            '#3498db', '#2ecc71', '#e74c3c', '#f39c12', '#9b59b6',
            '#1abc9c', '#34495e', '#e67e22', '#95a5a6', '#d35400'
        ];

        function showTab(tabId) {{
            document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            document.getElementById(tabId).classList.add('active');
            event.target.classList.add('active');
        }}

        // Runtime pie chart
        new Chart(document.getElementById('runtimeChart'), {{
            type: 'doughnut',
            data: {{
                labels: processNames,
                datasets: [{{
                    data: runtimeData,
                    backgroundColor: colors.concat(colors),
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Runtime Distribution (minutes)'
                    }}
                }}
            }}
        }});

        // Memory bar chart
        new Chart(document.getElementById('memoryChart'), {{
            type: 'bar',
            data: {{
                labels: processNames,
                datasets: [{{
                    label: 'Peak Memory (GB)',
                    data: memoryData,
                    backgroundColor: '#3498db',
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Peak Memory by Process'
                    }}
                }}
            }}
        }});

        // Task count chart
        new Chart(document.getElementById('taskCountChart'), {{
            type: 'bar',
            data: {{
                labels: processNames,
                datasets: [{{
                    label: 'Number of Tasks',
                    data: taskCounts,
                    backgroundColor: '#2ecc71',
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Task Count by Process'
                    }}
                }}
            }}
        }});

        // Memory detail chart
        new Chart(document.getElementById('memoryDetailChart'), {{
            type: 'bar',
            data: {{
                labels: processNames,
                datasets: [{{
                    label: 'Peak Memory (GB)',
                    data: memoryData,
                    backgroundColor: '#e74c3c',
                }}]
            }},
            options: {{
                indexAxis: 'y',
                responsive: true,
                maintainAspectRatio: false,
            }}
        }});

        // Runtime detail chart
        new Chart(document.getElementById('runtimeDetailChart'), {{
            type: 'bar',
            data: {{
                labels: processNames,
                datasets: [{{
                    label: 'Total Runtime (minutes)',
                    data: runtimeData,
                    backgroundColor: '#9b59b6',
                }}]
            }},
            options: {{
                indexAxis: 'y',
                responsive: true,
                maintainAspectRatio: false,
            }}
        }});

        // Efficiency chart
        new Chart(document.getElementById('efficiencyChart'), {{
            type: 'bar',
            data: {{
                labels: processNames,
                datasets: [{{
                    label: 'Memory Efficiency (%)',
                    data: efficiencyData,
                    backgroundColor: efficiencyData.map(e =>
                        e < 30 ? '#e74c3c' : e < 70 ? '#f1c40f' : '#2ecc71'
                    ),
                }}]
            }},
            options: {{
                indexAxis: 'y',
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    x: {{
                        max: 100,
                        title: {{
                            display: true,
                            text: 'Efficiency (Peak Used / Requested) %'
                        }}
                    }}
                }}
            }}
        }});
    </script>
</body>
</html>'''

    with open(output_path, 'w') as f:
        f.write(html)

    print(f"Report generated: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate resource usage visualization from Nextflow trace file'
    )
    parser.add_argument('--trace', required=True, help='Path to Nextflow trace TSV file')
    parser.add_argument('--outdir', default='.', help='Output directory for HTML report')
    parser.add_argument('--output', default='resource_usage.html', help='Output filename')

    args = parser.parse_args()

    trace_path = Path(args.trace)
    if not trace_path.exists():
        print(f"Error: Trace file not found: {trace_path}")
        sys.exit(1)

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / args.output

    print(f"Parsing trace file: {trace_path}")
    tasks = parse_trace_file(trace_path)
    print(f"Found {len(tasks)} tasks")

    processes = aggregate_by_process(tasks)
    print(f"Aggregated into {len(processes)} processes")

    generate_html_report(processes, output_path, str(trace_path))


if __name__ == '__main__':
    main()
