"""Export subgraph_pipeline SPN simulation to Obsidian vault."""

from __future__ import annotations

import csv
import json
import re
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from gin_tool.subgraph_pipeline.spn_policy import DEFAULT_SPN_POLICY, SPN_PIPELINE_VERSION

# inline dashboard (from scBionet obsidian_bridge)
_DASHBOARD_HTML_SRC = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>SPN dashboard</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
  <style>
    body { font-family: system-ui, sans-serif; margin: 1rem; background: #fafafa; }
    .row { display: flex; flex-wrap: wrap; gap: 1rem; }
    .panel { background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 1rem; flex: 1 1 420px; }
    select { margin: 0.25rem 0.5rem 0.25rem 0; }
    canvas { max-height: 360px; }
    #meta { font-size: 0.9rem; color: #444; }
  </style>
</head>
<body>
  <h1>SPN trajectory</h1>
  <p id="meta">loading…</p>
  <div class="row">
    <div class="panel">
      <h2>Place time series</h2>
      <label>Place <select id="placeSelect"></select></label>
      <canvas id="lineChart"></canvas>
    </div>
    <div class="panel">
      <h2>Top variable places (σ)</h2>
      <canvas id="barChart"></canvas>
    </div>
  </div>
  <script>
    const RUN = "__RUN_ID__";
    const base = `_meta/spn/runs/${RUN}/`;
    let traj = null, lineChart = null, barChart = null;
    function std(arr) {
      const m = arr.reduce((a,b)=>a+b,0)/arr.length;
      return Math.sqrt(arr.reduce((a,b)=>a+(b-m)**2,0)/arr.length);
    }
    async function load() {
      const resp = await fetch(base + 'trajectory.json');
      traj = await resp.json();
      const places = traj.places;
      const means = traj.mean_markings;
      document.getElementById('meta').textContent =
        `${places.length} places · ${means.length-1} blocks · ${traj.n_runs} runs`;
      const stats = places.map((p,i) => ({
        place: p, idx: i,
        series: means.map(row => row[i]),
        sigma: std(means.map(row => row[i]))
      })).sort((a,b)=>b.sigma-a.sigma);
      const sel = document.getElementById('placeSelect');
      stats.forEach(s => {
        const o = document.createElement('option');
        o.value = s.idx; o.textContent = s.place;
        sel.appendChild(o);
      });
      sel.onchange = () => drawLine(+sel.value, stats);
      drawLine(stats[0].idx, stats);
      drawBar(stats.slice(0, 20));
    }
    function drawLine(idx, stats) {
      const s = stats.find(x => x.idx === idx) || stats[0];
      const ctx = document.getElementById('lineChart');
      if (lineChart) lineChart.destroy();
      lineChart = new Chart(ctx, {
        type: 'line',
        data: { labels: s.series.map((_,i)=>i), datasets: [{ label: s.place, data: s.series, borderColor:'#2563eb', tension:0.15 }] },
        options: { scales: { y: { beginAtZero: true } } }
      });
    }
    function drawBar(top) {
      const ctx = document.getElementById('barChart');
      if (barChart) barChart.destroy();
      barChart = new Chart(ctx, {
        type: 'bar',
        data: {
          labels: top.map(t => t.place.length > 18 ? t.place.slice(0,16)+'…' : t.place),
          datasets: [{ label: 'σ(block mean)', data: top.map(t => t.sigma), backgroundColor:'#059669' }]
        },
        options: { indexAxis: 'y', plugins: { legend: { display: false } } }
      });
    }
    load().catch(e => { document.getElementById('meta').textContent = 'Error: ' + e; });
  </script>
</body>
</html>
"""


def safe_note_name(name: str) -> str:
    s = str(name).strip()
    s = s.replace("/", "_").replace("\\", "_").replace(":", "_")
    s = re.sub(r"[\x00-\x1f]", "", s)
    s = s.replace("|", "_").replace("#", "_").replace("^", "_").replace("[", "(").replace("]", ")")
    if not s:
        s = "unnamed"
    return s.rstrip(".")[:180]


def _topology_from_model(model) -> dict:
    return {
        "places": list(model.places),
        "interactions": [
            {
                "source": e.source,
                "target": e.target,
                "kind": e.kind,
                "transition_id": getattr(e, "transition_id", ""),
                "weight": getattr(e, "weight", 1.0),
            }
            for e in model.interactions
        ],
        "n_places": len(model.places),
        "n_interactions": len(model.interactions),
    }


def export_subgraph_spn_obsidian(
    subgraph_dir: str | Path,
    out_dir: str | Path,
    *,
    stem: str = "p53",
    write_place_notes: bool = True,
    top_neighbors: int = 15,
) -> dict:
    """Build Obsidian vault from subgraph_pipeline spn_sim output."""
    subgraph_dir = Path(subgraph_dir)
    out_dir = Path(out_dir)
    sim_dir = subgraph_dir / "spn_sim"
    traj_path = sim_dir / f"{stem}_trajectory.json"
    csv_path = sim_dir / f"{stem}_trajectory.csv"
    summary_path = sim_dir / f"{stem}_sim_summary.json"
    graphml_src = subgraph_dir / "spn.graphml"

    if not traj_path.is_file():
        raise FileNotFoundError(traj_path)

    summary = json.loads(summary_path.read_text(encoding="utf-8")) if summary_path.is_file() else {}
    traj = json.loads(traj_path.read_text(encoding="utf-8"))
    run_id = f"{stem}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"

    # topology from live model rebuild
    from gin_tool.subgraph_pipeline.spn_simulate import (
        build_petri_net_from_bipartite_edges,
        load_gin_edges_tsv,
        petri_net_to_spn_model,
        seed_gin_names_from_manifest,
    )
    edges_path = subgraph_dir / "gin_edges.tsv"
    if not edges_path.is_file():
        edges_path = subgraph_dir / "subgraph_edges.tsv"
    edges = load_gin_edges_tsv(edges_path)
    manifest = subgraph_dir / "seed_manifest.json"
    seeds = seed_gin_names_from_manifest(manifest) if manifest.is_file() else set()
    net = build_petri_net_from_bipartite_edges(edges, seed_gin_names=seeds)
    model = petri_net_to_spn_model(net)
    topo = _topology_from_model(model)

    meta_run = out_dir / "_meta" / "spn" / "runs" / run_id
    meta_root = out_dir / "_meta" / "spn"
    nodes_dir = out_dir / "nodes"
    meta_run.mkdir(parents=True, exist_ok=True)
    meta_root.mkdir(parents=True, exist_ok=True)
    if write_place_notes and nodes_dir.exists():
        shutil.rmtree(nodes_dir)
    if write_place_notes:
        nodes_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(traj_path, meta_run / "trajectory.json")
    if csv_path.is_file():
        shutil.copy2(csv_path, meta_run / "trajectory.csv")
    (meta_run / "topology.json").write_text(json.dumps(topo, indent=2), encoding="utf-8")
    if graphml_src.is_file():
        shutil.copy2(graphml_src, meta_run / "spn.graphml")
        shutil.copy2(graphml_src, out_dir / "spn.graphml")

    manifest_out = {
        "run_id": run_id,
        "tag": stem,
        "spn_method_version": summary.get("spn_pipeline_version")
        or summary.get("firing_policy", DEFAULT_SPN_POLICY.version),
        "graphml": str(graphml_src),
        "simulator": {"blocks": summary.get("blocks"), "runs": summary.get("runs")},
        "topology": {"n_places": topo["n_places"], "n_interactions": topo["n_interactions"]},
        "outputs": {
            "trajectory_json": str(meta_run / "trajectory.json"),
            "topology_json": str(meta_run / "topology.json"),
        },
    }
    (meta_run / "run_manifest.json").write_text(json.dumps(manifest_out, indent=2), encoding="utf-8")

    neighbors: dict[str, set[str]] = defaultdict(set)
    for edge in topo["interactions"]:
        src, tgt = edge["source"], edge["target"]
        if src != tgt:
            neighbors[src].add(tgt)
            neighbors[tgt].add(src)

    places: list[str] = list(traj["places"])
    means = traj["mean_markings"]
    n_blocks = len(means)
    summary_rows: list[dict] = []

    for i, place in enumerate(places):
        series = [float(means[b][i]) for b in range(n_blocks)]
        note_id = safe_note_name(place)
        start, end = series[0], series[-1]
        mx = max(series)
        sigma = (sum((x - sum(series) / len(series)) ** 2 for x in series) / len(series)) ** 0.5
        summary_rows.append(
            {
                "id": note_id,
                "place": place,
                "mean_start": f"{start:.4f}",
                "mean_end": f"{end:.4f}",
                "delta": f"{end - start:.4f}",
                "max": f"{mx:.4f}",
                "sigma": f"{sigma:.4f}",
            }
        )
        if not write_place_notes:
            continue
        nbs = sorted(neighbors.get(place, []), key=safe_note_name)[:top_neighbors]
        lines = [
            "---",
            "tags: [spn-place]",
            "spn_place: true",
            f'run_id: "{run_id}"',
            f"mean_start: {start:.4f}",
            f"mean_end: {end:.4f}",
            f"delta: {end - start:.4f}",
            f"max_token: {mx:.4f}",
            f"sigma: {sigma:.4f}",
            f"blocks: {n_blocks - 1}",
            "---",
            "",
            f"# {place}",
            "",
            f"- run: `{run_id}`",
            f"- mean token: {start:.3g} → {end:.3g} (Δ {end - start:+.3g})",
            f"- max: {mx:.3g} · σ: {sigma:.3g}",
            "",
        ]
        if nbs:
            lines += ["## Neighbors (SPN topology)", ""]
            for nb in nbs:
                lines.append(f"- [[{safe_note_name(nb)}|{nb}]]")
            lines.append("")
        (nodes_dir / f"{note_id}.md").write_text("\n".join(lines), encoding="utf-8")

    summary_rows.sort(key=lambda r: -float(r["sigma"]))
    with (meta_root / "places_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "place", "mean_start", "mean_end", "delta", "max", "sigma"])
        w.writeheader()
        w.writerows(summary_rows)

    config = {
        "run_id": run_id,
        "tag": stem,
        "spn_method_version": summary.get("spn_pipeline_version")
        or summary.get("firing_policy", SPN_PIPELINE_VERSION),
        "graphml": "spn.graphml",
        "vault_paths": {
            "trajectory_json": f"_meta/spn/runs/{run_id}/trajectory.json",
            "topology_json": f"_meta/spn/runs/{run_id}/topology.json",
            "graphml": f"_meta/spn/runs/{run_id}/spn.graphml",
        },
    }
    (meta_root / "config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    html = _DASHBOARD_HTML_SRC.replace("__RUN_ID__", run_id)
    (out_dir / "spn_dashboard.html").write_text(html, encoding="utf-8")

    dashboard_md = f"""---
tags: [spn-dashboard]
run_id: "{run_id}"
---

# SPN dashboard · {stem}

Open [[spn_dashboard.html]] in browser (Chart.js reads `_meta/spn/runs/{run_id}/trajectory.json`).

GraphML: [[spn.graphml]] · BioLayout Express3D

## Dataview

```dataview
TABLE mean_start, mean_end, delta, max_token, sigma
FROM "nodes"
WHERE spn_place = true
SORT sigma DESC
LIMIT 30
```

## Run

- firing: `{summary.get("firing_policy")}`
- places: {topo["n_places"]} · interactions: {topo["n_interactions"]}
- blocks: {summary.get("blocks")} · runs: {summary.get("runs")}
- top token: {summary.get("top_places_final_block", [[None, 0]])[0][1] if summary.get("top_places_final_block") else "n/a"}
"""
    (out_dir / "dashboard.md").write_text(dashboard_md, encoding="utf-8")

    export_summary = {
        "run_id": run_id,
        "out_dir": str(out_dir.resolve()),
        "n_place_notes": len(summary_rows) if write_place_notes else 0,
        "graphml": str((out_dir / "spn.graphml").resolve()) if graphml_src.is_file() else None,
        "dashboard_md": str((out_dir / "dashboard.md").resolve()),
        "dashboard_html": str((out_dir / "spn_dashboard.html").resolve()),
    }
    (meta_root / "export_summary.json").write_text(json.dumps(export_summary, indent=2), encoding="utf-8")
    return export_summary
