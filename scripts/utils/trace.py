#!/usr/bin/env python3
"""Trace file aggregation and resource usage visualization utilities."""

import re
from datetime import datetime

from .parsers import parse_duration, parse_memory, format_bytes, format_duration_str


def parse_timestamp(ts_str):
    """Parse Nextflow timestamp string to datetime."""
    if not ts_str or ts_str == '-':
        return None
    try:
        # Nextflow format: 2024-01-12 23:45:30.123
        return datetime.strptime(ts_str.split('.')[0], '%Y-%m-%d %H:%M:%S')
    except (ValueError, AttributeError):
        return None


def aggregate_trace_by_process(tasks):
    """Aggregate task statistics by process name."""
    processes = {}
    for task in tasks:
        name = task.get('name', 'unknown')
        # Remove parameters in parentheses, handling nested parens like "(Streptomyces coelicolor A3(2))"
        base_name = re.sub(r'\s*\([^()]*(?:\([^()]*\)[^()]*)*\)$', '', name).strip()
        if base_name not in processes:
            processes[base_name] = {
                'name': base_name, 'count': 0, 'completed': 0, 'failed': 0,
                'total_runtime': 0, 'peak_memory': 0, 'total_memory': 0,
                'peak_cpu': 0, 'total_cpu': 0, 'memory_requested': 0,
            }
        proc = processes[base_name]
        proc['count'] += 1
        status = task.get('status', '')
        if status in ('COMPLETED', 'CACHED'):
            proc['completed'] += 1
        elif status in ('FAILED', 'ABORTED'):
            proc['failed'] += 1
        runtime = parse_duration(task.get('realtime', '0'))
        peak_rss = parse_memory(task.get('peak_rss', '0'))
        cpu_str = task.get('%cpu', '0').replace('%', '')
        cpu_pct = float(cpu_str) if cpu_str and cpu_str != '-' else 0
        mem_requested = parse_memory(task.get('memory', '0'))
        proc['total_runtime'] += runtime
        proc['total_memory'] += peak_rss
        proc['total_cpu'] += cpu_pct
        proc['peak_memory'] = max(proc['peak_memory'], peak_rss)
        proc['peak_cpu'] = max(proc['peak_cpu'], cpu_pct)
        proc['memory_requested'] = max(proc['memory_requested'], mem_requested)

    for proc in processes.values():
        if proc['count'] > 0:
            proc['avg_runtime'] = proc['total_runtime'] / proc['count']
            proc['avg_memory'] = proc['total_memory'] / proc['count']
            proc['avg_cpu'] = proc['total_cpu'] / proc['count']
    return processes


def generate_gantt_chart_html(tasks):
    """Generate an HTML/CSS-based Gantt chart showing task execution timeline.

    Aggregates tasks by module/process name to show one bar per module instead
    of individual genome tasks.
    """
    if not tasks:
        return ''

    # Parse timestamps and filter tasks with valid timing data
    timed_tasks = []
    for task in tasks:
        start = parse_timestamp(task.get('start', ''))
        complete = parse_timestamp(task.get('complete', ''))
        if start and complete:
            # Extract base process name (remove genome-specific parts in parentheses)
            full_name = task.get('name', 'Unknown')
            base_name = re.sub(r'\s*\([^()]*(?:\([^()]*\)[^()]*)*\)$', '', full_name).strip()
            timed_tasks.append({
                'name': full_name,
                'base_name': base_name,
                'status': task.get('status', ''),
                'start': start,
                'complete': complete,
                'duration': (complete - start).total_seconds()
            })

    if not timed_tasks:
        return '<p style="color: #7f8c8d; font-style: italic;">Timeline data not available. Run pipeline without -resume to generate timing data.</p>'

    # Aggregate tasks by module (base process name)
    modules = {}
    for task in timed_tasks:
        base_name = task['base_name']
        if base_name not in modules:
            modules[base_name] = {
                'name': base_name,
                'tasks': [],
                'start': task['start'],
                'complete': task['complete'],
                'total_duration': 0,
                'failed': 0,
                'succeeded': 0
            }
        mod = modules[base_name]
        mod['tasks'].append(task)
        mod['start'] = min(mod['start'], task['start'])
        mod['complete'] = max(mod['complete'], task['complete'])
        mod['total_duration'] += task['duration']
        if task['status'] in ('FAILED', 'ABORTED'):
            mod['failed'] += 1
        else:
            mod['succeeded'] += 1

    # Convert to list and sort by start time
    module_list = sorted(modules.values(), key=lambda x: x['start'])

    # Calculate timeline bounds
    min_time = min(m['start'] for m in module_list)
    max_time = max(m['complete'] for m in module_list)
    total_duration = (max_time - min_time).total_seconds()

    if total_duration <= 0:
        return '<p style="color: #7f8c8d; font-style: italic;">Timeline duration too short to display.</p>'

    # Color palette for different process types
    color_palette = ['#3498db', '#2ecc71', '#9b59b6', '#f39c12', '#1abc9c', '#e67e22', '#34495e', '#16a085']
    process_colors = {}
    for i, mod in enumerate(module_list):
        process_colors[mod['name']] = color_palette[i % len(color_palette)]

    # Build Gantt chart HTML
    html = '''
    <div class="gantt-container" style="margin: 20px 0; overflow-x: auto;">
        <div class="gantt-chart" style="position: relative; min-width: 600px; background: #f8f9fa; border-radius: 8px; padding: 15px;">
            <div class="gantt-header" style="display: flex; justify-content: space-between; margin-bottom: 10px; font-size: 0.85em; color: #666;">
                <span>Start: ''' + min_time.strftime('%H:%M:%S') + '''</span>
                <span>Total Duration: ''' + format_duration_str(total_duration) + '''</span>
                <span>End: ''' + max_time.strftime('%H:%M:%S') + '''</span>
            </div>
            <div class="gantt-timeline" style="position: relative; border-left: 2px solid #ddd; border-right: 2px solid #ddd;">
    '''

    # Render each module as a bar
    for mod in module_list:
        offset_seconds = (mod['start'] - min_time).total_seconds()
        span_duration = (mod['complete'] - mod['start']).total_seconds()
        left_pct = (offset_seconds / total_duration) * 100
        width_pct = (span_duration / total_duration) * 100
        width_pct = max(width_pct, 1.0)  # Minimum width for visibility

        color = process_colors.get(mod['name'], '#3498db')
        task_count = len(mod['tasks'])

        # Show failure indicator if any tasks failed
        if mod['failed'] > 0:
            border = '2px solid #e74c3c'
            status_text = f"{mod['succeeded']} succeeded, {mod['failed']} failed"
        else:
            border = 'none'
            status_text = f"{task_count} task{'s' if task_count > 1 else ''}"

        short_name = mod['name'].split(':')[-1] if ':' in mod['name'] else mod['name']
        span_str = format_duration_str(span_duration)
        cpu_str = format_duration_str(mod['total_duration'])

        html += f'''
                <div class="gantt-bar" style="position: relative; height: 32px; margin: 6px 0;">
                    <div style="position: absolute; left: {left_pct:.2f}%; width: {width_pct:.2f}%; height: 100%;
                                background: {color}; border-radius: 4px; border: {border};
                                display: flex; align-items: center; padding: 0 10px; overflow: hidden;
                                box-shadow: 0 1px 3px rgba(0,0,0,0.15);"
                         title="{mod['name']}&#10;{status_text}&#10;Wall time: {span_str}&#10;CPU time: {cpu_str}">
                        <span style="color: white; font-size: 0.8em; font-weight: 500; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
                            {short_name} ({task_count})
                        </span>
                    </div>
                </div>
        '''

    html += '''
            </div>
        </div>
    </div>

    <!-- Legend -->
    <div style="display: flex; flex-wrap: wrap; gap: 15px; margin-top: 10px; font-size: 0.85em;">
    '''

    for mod in module_list:
        color = process_colors[mod['name']]
        short_name = mod['name'].split(':')[-1] if ':' in mod['name'] else mod['name']
        task_count = len(mod['tasks'])
        html += f'''
        <div style="display: flex; align-items: center; gap: 5px;">
            <div style="width: 16px; height: 16px; background: {color}; border-radius: 3px;"></div>
            <span>{short_name} ({task_count})</span>
        </div>
        '''

    html += '</div>'

    return html


def generate_resource_usage_html(tasks, processes):
    """Generate HTML content for resource usage tab."""
    if not tasks:
        return '''
        <div class="info-box warning">
            <p>No resource usage data available.</p>
            <p>Resource tracking data will appear here after running the pipeline.</p>
        </div>
        '''

    sorted_procs = sorted(processes.values(), key=lambda x: x['total_runtime'], reverse=True)
    total_tasks = sum(p['count'] for p in sorted_procs)
    total_completed = sum(p['completed'] for p in sorted_procs)
    total_failed = sum(p['failed'] for p in sorted_procs)
    total_runtime = sum(p['total_runtime'] for p in sorted_procs)
    peak_memory_overall = max((p['peak_memory'] for p in sorted_procs), default=0)

    # Generate Gantt chart
    gantt_html = generate_gantt_chart_html(tasks)

    html = f'''
    <div class="stats-container" style="margin-bottom: 25px;">
        <div class="stat-box"><div class="stat-value">{total_tasks}</div><div class="stat-label">Total Tasks</div></div>
        <div class="stat-box highlight"><div class="stat-value">{total_completed}</div><div class="stat-label">Completed</div></div>
        <div class="stat-box{"" if total_failed == 0 else " warning"}"><div class="stat-value">{total_failed}</div><div class="stat-label">Failed</div></div>
        <div class="stat-box"><div class="stat-value">{format_duration_str(total_runtime)}</div><div class="stat-label">Total Runtime</div></div>
        <div class="stat-box"><div class="stat-value">{format_bytes(peak_memory_overall)}</div><div class="stat-label">Peak Memory</div></div>
    </div>

    <h3>Execution Timeline</h3>
    {gantt_html}

    <h3>Resource Usage by Process</h3>
    <table class="data-table" style="margin-bottom: 25px;">
        <thead>
            <tr>
                <th>Process</th>
                <th>Tasks</th>
                <th>Completed</th>
                <th>Failed</th>
                <th>Total Runtime</th>
                <th>Avg Runtime</th>
                <th>Peak Memory</th>
                <th>Avg CPU</th>
            </tr>
        </thead>
        <tbody>
    '''

    for p in sorted_procs:
        fail_style = ' style="color: #c0392b; font-weight: bold;"' if p['failed'] > 0 else ''
        html += f'''
            <tr>
                <td><strong>{p['name']}</strong></td>
                <td>{p['count']}</td>
                <td style="color: #27ae60;">{p['completed']}</td>
                <td{fail_style}>{p['failed']}</td>
                <td>{format_duration_str(p['total_runtime'])}</td>
                <td>{format_duration_str(p.get('avg_runtime', 0))}</td>
                <td>{format_bytes(p['peak_memory'])}</td>
                <td>{p.get('avg_cpu', 0):.1f}%</td>
            </tr>
        '''

    html += '''
        </tbody>
    </table>
    '''

    return html
