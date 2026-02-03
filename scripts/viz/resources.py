#!/usr/bin/env python3
"""Resource usage visualization functions."""

import re
from datetime import datetime
import sys
sys.path.insert(0, str(__file__).rsplit('/', 2)[0])
from utils.parsers import parse_duration, parse_memory, format_bytes, format_duration_str


def aggregate_trace_by_process(tasks):
    """Aggregate task statistics by process name."""
    processes = {}
    for task in tasks:
        name = task.get('name', 'unknown')
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


def parse_timestamp(ts_str):
    """Parse Nextflow timestamp string to datetime."""
    if not ts_str or ts_str == '-':
        return None
    try:
        return datetime.strptime(ts_str.split('.')[0], '%Y-%m-%d %H:%M:%S')
    except ValueError:
        return None


def generate_gantt_chart_html(tasks):
    """Generate Gantt chart HTML for task execution timeline."""
    valid_tasks = []
    for task in tasks:
        start = parse_timestamp(task.get('start'))
        complete = parse_timestamp(task.get('complete'))
        if start and complete:
            valid_tasks.append({
                'name': task.get('name', 'Unknown'),
                'start': start,
                'end': complete,
                'status': task.get('status', '')
            })

    if not valid_tasks:
        return '<div class="info-box"><p>No timeline data available.</p></div>'

    valid_tasks.sort(key=lambda x: x['start'])
    min_time = min(t['start'] for t in valid_tasks)
    max_time = max(t['end'] for t in valid_tasks)
    total_seconds = (max_time - min_time).total_seconds()

    if total_seconds <= 0:
        return '<div class="info-box"><p>Timeline too short to display.</p></div>'

    # Group by process
    process_tasks = {}
    for task in valid_tasks:
        name = re.sub(r'\s*\([^()]*(?:\([^()]*\)[^()]*)*\)$', '', task['name']).strip()
        if name not in process_tasks:
            process_tasks[name] = []
        process_tasks[name].append(task)

    html = '''
    <div style="position: relative; width: 100%; height: auto; margin-bottom: 20px;">
        <div style="display: flex; flex-direction: column; gap: 4px;">
    '''

    colors = {'COMPLETED': '#27ae60', 'CACHED': '#3498db', 'FAILED': '#e74c3c', 'ABORTED': '#e74c3c'}

    for proc_name, proc_tasks in sorted(process_tasks.items()):
        html += f'''
        <div style="display: flex; align-items: center; height: 24px;">
            <div style="width: 200px; font-size: 11px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; padding-right: 10px;" title="{proc_name}">{proc_name}</div>
            <div style="flex: 1; position: relative; height: 100%; background: #f0f0f0; border-radius: 2px;">
        '''
        for task in proc_tasks:
            left = ((task['start'] - min_time).total_seconds() / total_seconds) * 100
            width = ((task['end'] - task['start']).total_seconds() / total_seconds) * 100
            width = max(width, 0.5)
            color = colors.get(task['status'], '#95a5a6')
            html += f'''
                <div style="position: absolute; left: {left:.2f}%; width: {width:.2f}%; height: 100%; background: {color}; border-radius: 2px;" title="{task['name']}"></div>
            '''
        html += '''
            </div>
        </div>
        '''

    html += '''
        </div>
        <div style="display: flex; justify-content: space-between; margin-top: 8px; font-size: 10px; color: #666;">
    '''
    html += f'<span>{min_time.strftime("%H:%M:%S")}</span>'
    html += f'<span>{max_time.strftime("%H:%M:%S")}</span>'
    html += '</div></div>'

    # Legend
    html += '''
    <div style="margin-top: 10px; font-size: 11px;">
        <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #27ae60; border-radius: 2px; vertical-align: middle; margin-right: 4px;"></span>Completed</span>
        <span style="margin-right: 15px;"><span style="display: inline-block; width: 12px; height: 12px; background: #3498db; border-radius: 2px; vertical-align: middle; margin-right: 4px;"></span>Cached</span>
        <span><span style="display: inline-block; width: 12px; height: 12px; background: #e74c3c; border-radius: 2px; vertical-align: middle; margin-right: 4px;"></span>Failed</span>
    </div>
    '''

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
