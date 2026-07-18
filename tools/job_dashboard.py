# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "fastapi",
#     "uvicorn",
# ]
# ///

import argparse
import json
import os
import sys
from pathlib import Path

import uvicorn
from fastapi import FastAPI
from fastapi.responses import HTMLResponse

app = FastAPI()
target_dir = ""
cached_total_target_cases = None


def get_stats():
    global cached_total_target_cases
    stats = {
        "total": 0,
        "success": 0,
        "failed": 0,
        "pending": 0,
        "remaining": 0,
        "accuracy": 0.0,
        "percent_complete": 0.0,
        "jobs": [],
    }

    registry_path = Path(target_dir) / "case_registry.json"
    summary_path = Path(target_dir) / "summary_roundtrip.json"

    if cached_total_target_cases is None:
        total_target_cases = 0
        if registry_path.exists():
            try:
                with open(registry_path, "r") as f:
                    registry_data = json.load(f)
                    total_target_cases = len(registry_data)
            except Exception:
                pass
        cached_total_target_cases = total_target_cases
    else:
        total_target_cases = cached_total_target_cases

    all_jobs = []

    try:
        # summary_roundtrip.json is the live source the sweep appends to on every
        # molecule, so it is authoritative for "processed so far". case_registry.json
        # is a derived classification that lags behind it (and may be absent), so it
        # is only a fallback. Either way, "total" is the count of processed cases.
        data = None
        if summary_path.exists():
            try:
                with open(summary_path, "r") as f:
                    data = json.load(f)
                label_key = "status"
            except json.JSONDecodeError:
                try:
                    with open(summary_path, "r") as f:
                        raw = f.read()
                    last_brace = raw.rfind("}")
                    if last_brace != -1:
                        data = json.loads(raw[: last_brace + 1] + "]")
                        label_key = "status"
                except Exception:
                    pass

        if data is None and registry_path.exists():
            try:
                with open(registry_path, "r") as f:
                    data = json.load(f)
                label_key = "class"
            except json.JSONDecodeError:
                pass

        if data is None:
            data = []
            label_key = "status"

        # Deduplicate by molecule name, keeping the latest attempt
        unique_jobs = {}
        for item in data:
            mol = item.get("molecule")
            if mol:
                unique_jobs[mol] = item
            else:
                unique_jobs[id(item)] = item
        data = list(unique_jobs.values())

        stats["total"] = len(data)
        stats["remaining"] = max(0, total_target_cases - stats["total"])
        for item in data:
            # A case is pending only while still in flight (label begins with
            # "pending", e.g. "pending_g-xtb"). Every other non-success label --
            # "failed", "error", and the registry's granular classes like
            # "no_conformers" or "high_rmsd" -- is a real failure, not pending.
            label = item.get(label_key, "pending")
            if label == "success":
                stats["success"] += 1
            elif label.startswith("pending"):
                stats["pending"] += 1
            else:
                stats["failed"] += 1
        if stats["total"] > 0:
            stats["accuracy"] = (stats["success"] / stats["total"]) * 100
        else:
            stats["accuracy"] = 0.0

        if total_target_cases > 0:
            stats["percent_complete"] = min(100.0, (stats["total"] / total_target_cases) * 100)
        else:
            stats["percent_complete"] = 100.0 if stats["total"] > 0 else 0.0

        all_jobs = data

        # Add a derived timestamp to each job for robust sorting
        import datetime

        for job in all_jobs:
            if job.get("saved_at"):
                job["_sort_time"] = job["saved_at"]
            else:
                # Try to get file mod time if available
                input_file = job.get("input_xyz")
                if input_file and os.path.exists(input_file):
                    mtime = os.path.getmtime(input_file)
                    job["_sort_time"] = datetime.datetime.fromtimestamp(mtime).isoformat()
                else:
                    job["_sort_time"] = "1970-01-01T00:00:00"

        # Sort explicitly by timestamp, newest first
        all_jobs.sort(key=lambda x: x.get("_sort_time", ""), reverse=True)
        stats["jobs"] = all_jobs[:50]

        # Clean up internal sort key before sending to frontend
        for job in stats["jobs"]:
            if "_sort_time" in job:
                job["saved_at"] = job.pop("_sort_time")  # Use as display timestamp

    except Exception as e:
        print(f"Error reading stats: {e}")

    return stats


@app.get("/api/data")
async def data_endpoint():
    return get_stats()


@app.get("/", response_class=HTMLResponse)
async def index():
    return """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>OIN-SMILES Roundtrip Sweep Dashboard</title>
        <link rel="preconnect" href="https://fonts.googleapis.com">
        <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
        <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap" rel="stylesheet">
        <style>
            :root {
                --bg-gradient-start: #0f172a;
                --bg-gradient-end: #1e1b4b;
                --glass-bg: rgba(255, 255, 255, 0.05);
                --glass-border: rgba(255, 255, 255, 0.1);
                --text-main: #f8fafc;
                --text-muted: #94a3b8;
                --success: #10b981;
                --failed: #ef4444;
                --pending: #f59e0b;
                --card-radius: 16px;
            }

            * {
                box-sizing: border-box;
                margin: 0;
                padding: 0;
            }

            body {
                font-family: 'Inter', sans-serif;
                background: linear-gradient(135deg, var(--bg-gradient-start), var(--bg-gradient-end));
                color: var(--text-main);
                min-height: 100vh;
                padding: 2rem;
                display: flex;
                flex-direction: column;
                align-items: center;
                background-attachment: fixed;
            }

            .container {
                width: 100%;
                max-width: 1200px;
                display: flex;
                flex-direction: column;
                gap: 2rem;
            }

            header {
                display: flex;
                justify-content: space-between;
                align-items: center;
                padding-bottom: 1rem;
                border-bottom: 1px solid var(--glass-border);
            }

            header h1 {
                font-weight: 700;
                font-size: 2.5rem;
                background: -webkit-linear-gradient(#fff, #cbd5e1);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.02em;
            }

            .status-indicator {
                display: flex;
                align-items: center;
                gap: 0.5rem;
                font-size: 0.875rem;
                color: var(--success);
                background: var(--glass-bg);
                padding: 0.5rem 1rem;
                border-radius: 9999px;
                border: 1px solid rgba(16, 185, 129, 0.2);
            }

            .status-dot {
                width: 8px;
                height: 8px;
                background-color: var(--success);
                border-radius: 50%;
                box-shadow: 0 0 8px var(--success);
                animation: pulse 2s infinite;
            }

            @keyframes pulse {
                0% { box-shadow: 0 0 0 0 rgba(16, 185, 129, 0.7); }
                70% { box-shadow: 0 0 0 6px rgba(16, 185, 129, 0); }
                100% { box-shadow: 0 0 0 0 rgba(16, 185, 129, 0); }
            }

            .stats-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 1.5rem;
            }

            .stat-card {
                background: var(--glass-bg);
                backdrop-filter: blur(12px);
                -webkit-backdrop-filter: blur(12px);
                border: 1px solid var(--glass-border);
                border-radius: var(--card-radius);
                padding: 1.5rem;
                display: flex;
                flex-direction: column;
                gap: 0.5rem;
                transition: transform 0.3s ease, box-shadow 0.3s ease;
            }

            .stat-card:hover {
                transform: translateY(-5px);
                box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.5);
            }

            .stat-title {
                font-size: 0.875rem;
                text-transform: uppercase;
                letter-spacing: 0.05em;
                color: var(--text-muted);
                font-weight: 600;
            }

            .stat-value {
                font-size: 3rem;
                font-weight: 700;
                line-height: 1;
                white-space: nowrap;
            }

            .percent-complete {
                font-size: 1.25rem;
                font-weight: 500;
                color: var(--text-muted);
                vertical-align: bottom;
                margin-left: 4px;
            }

            .val-total { color: #fff; }
            .val-success { color: var(--success); text-shadow: 0 0 20px rgba(16,185,129,0.3); }
            .val-failed { color: var(--failed); text-shadow: 0 0 20px rgba(239,68,68,0.3); }
            .val-pending { color: var(--pending); text-shadow: 0 0 20px rgba(245,158,11,0.3); }
            .val-accuracy { color: #38bdf8; text-shadow: 0 0 20px rgba(56,189,248,0.3); }

            .jobs-section {
                background: var(--glass-bg);
                backdrop-filter: blur(12px);
                -webkit-backdrop-filter: blur(12px);
                border: 1px solid var(--glass-border);
                border-radius: var(--card-radius);
                padding: 1.5rem;
                overflow: hidden;
            }

            .jobs-section h2 {
                font-size: 1.5rem;
                margin-bottom: 1rem;
                color: var(--text-main);
            }

            .table-container {
                overflow-x: auto;
            }

            table {
                width: 100%;
                border-collapse: collapse;
                text-align: left;
            }

            th, td {
                padding: 1rem;
                border-bottom: 1px solid var(--glass-border);
            }

            th {
                color: var(--text-muted);
                font-weight: 600;
                font-size: 0.875rem;
                text-transform: uppercase;
                letter-spacing: 0.05em;
            }

            tr {
                transition: background-color 0.2s ease;
            }

            tr:hover {
                background-color: rgba(255, 255, 255, 0.03);
            }

            tr:last-child td {
                border-bottom: none;
            }

            .badge {
                padding: 0.25rem 0.75rem;
                border-radius: 9999px;
                font-size: 0.75rem;
                font-weight: 600;
                text-transform: uppercase;
            }

            .badge-success {
                background: rgba(16, 185, 129, 0.2);
                color: var(--success);
                border: 1px solid rgba(16, 185, 129, 0.3);
            }

            .badge-failed {
                background: rgba(239, 68, 68, 0.2);
                color: var(--failed);
                border: 1px solid rgba(239, 68, 68, 0.3);
            }

            .badge-pending {
                background: rgba(245, 158, 11, 0.2);
                color: var(--pending);
                border: 1px solid rgba(245, 158, 11, 0.3);
            }

            .fade-in {
                animation: fadeIn 0.5s ease-out forwards;
            }

            @keyframes fadeIn {
                from { opacity: 0; transform: translateY(10px); }
                to { opacity: 1; transform: translateY(0); }
            }

            .empty-state {
                text-align: center;
                padding: 3rem;
                color: var(--text-muted);
                font-style: italic;
            }

            .refresh-select {
                background: var(--glass-bg);
                color: var(--text-main);
                border: 1px solid var(--glass-border);
                padding: 0.5rem 1rem;
                border-radius: 9999px;
                outline: none;
                cursor: pointer;
                font-family: inherit;
                font-size: 0.875rem;
            }
            .refresh-select option {
                background: var(--bg-gradient-start);
                color: var(--text-main);
            }
        </style>
    </head>
    <body>
        <div class="container fade-in">
            <header>
                <h1>OIN-SMILES Roundtrip Sweep Dashboard</h1>
                <div style="display: flex; gap: 1rem; align-items: center;">
                    <select id="refresh-rate" class="refresh-select">
                        <option value="1000" selected>Refresh: 1 sec</option>
                        <option value="30000">Refresh: 30 sec</option>
                        <option value="60000">Refresh: 1 min</option>
                        <option value="300000">Refresh: 5 min</option>
                        <option value="600000">Refresh: 10 min</option>
                    </select>
                    <div class="status-indicator" id="connection-status">
                        <div class="status-dot"></div>
                        Live Updates Active
                    </div>
                </div>
            </header>

            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-title">Total Cases</div>
                    <div class="stat-value val-total" id="stat-total">0</div>
                </div>
                <div class="stat-card">
                    <div class="stat-title">Success</div>
                    <div class="stat-value val-success" id="stat-success">0</div>
                </div>
                <div class="stat-card">
                    <div class="stat-title">Failed / Errors</div>
                    <div class="stat-value val-failed" id="stat-failed">0</div>
                </div>
                <div class="stat-card">
                    <div class="stat-title">% Accuracy</div>
                    <div class="stat-value val-accuracy" id="stat-accuracy">0.0%</div>
                </div>
            </div>

            <div class="jobs-section">
                <h2>Recent Jobs</h2>
                <div class="table-container">
                    <table>
                        <thead>
                            <tr>
                                <th>Molecule</th>
                                <th>Status</th>
                                <th>Timestamp / Info</th>
                            </tr>
                        </thead>
                        <tbody id="jobs-tbody">
                            <tr>
                                <td colspan="3" class="empty-state">Loading data...</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <script>
            function updateDashboard() {
                fetch('/api/data')
                    .then(res => res.json())
                    .then(data => {
                        document.getElementById('stat-total').innerHTML = `${data.total} <span class="percent-complete">(${data.percent_complete.toFixed(1)}%)</span>`;
                        document.getElementById('stat-success').textContent = data.success;
                        document.getElementById('stat-failed').textContent = data.failed;

                        document.getElementById('stat-accuracy').textContent = data.accuracy.toFixed(1) + '%';

                        const tbody = document.getElementById('jobs-tbody');
                        tbody.innerHTML = '';

                        if (!data.jobs || data.jobs.length === 0) {
                            tbody.innerHTML = `<tr><td colspan="3" class="empty-state">No jobs found in directory.</td></tr>`;
                            return;
                        }

                        data.jobs.forEach(job => {
                            const tr = document.createElement('tr');

                            // Determine status
                            let status = job.class || job.status || 'pending';
                            let badgeClass = 'badge-pending';
                            let statusText = status.toUpperCase();

                            if (status === 'success') badgeClass = 'badge-success';
                            else if (status === 'error' || status === 'failed') badgeClass = 'badge-failed';

                            // Info column
                            let info = job.saved_at || job.tier_passed || job.commit_id || '-';

                            tr.innerHTML = `
                                <td style="font-weight: 600;">${job.molecule || 'Unknown'}</td>
                                <td><span class="badge ${badgeClass}">${statusText}</span></td>
                                <td style="color: var(--text-muted); font-size: 0.9em;">${info}</td>
                            `;
                            tbody.appendChild(tr);
                        });
                    })
                    .catch(err => {
                        console.error("Error fetching data:", err);
                        document.getElementById('connection-status').innerHTML = '<div style="width:8px;height:8px;background:red;border-radius:50%;"></div> Connection Error';
                        document.getElementById('connection-status').style.color = '#ef4444';
                        document.getElementById('connection-status').style.borderColor = 'rgba(239, 68, 68, 0.2)';
                    });
            }

            // Initial fetch
            updateDashboard();

            // Setup Polling
            let pollInterval = setInterval(updateDashboard, 1000);
            document.getElementById('refresh-rate').addEventListener('change', function(e) {
                clearInterval(pollInterval);
                const rate = parseInt(e.target.value, 10);
                pollInterval = setInterval(updateDashboard, rate);
            });
        </script>
    </body>
    </html>
    """


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Job Dashboard")
    parser.add_argument("folder", type=str, help="Folder to monitor")
    args = parser.parse_args()

    target_dir = args.folder
    if not os.path.isdir(target_dir):
        print(f"Error: {target_dir} is not a directory.")
        sys.exit(1)

    print(f"Monitoring {target_dir} for changes...")
    print("Starting server at http://localhost:8000")
    uvicorn.run(app, host="0.0.0.0", port=8000)
