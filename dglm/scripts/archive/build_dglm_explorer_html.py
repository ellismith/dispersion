#!/usr/bin/env python3
"""
Assembles the final dglm_explorer.html from the static template (head, CSS,
tab structure, and JS logic -- reproduced from the last known-working build)
plus a freshly-generated DATA JSON blob (built by build_dglm_explorer_data.R).

Model dropdown is collapsed to a single option (v3_sqrt_louvain) -- the new
pipeline only ever produces sqrt-shat, louvain-resolution output, so the old
4-way raw/sqrt x celltype/louvain comparison no longer applies.

Usage:
    python build_dglm_explorer_html.py --data /path/to/dglm_explorer_data.json --out /path/to/dglm_explorer.html
"""
import argparse
import json

TEMPLATE_HEAD = r'''<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>DGLM Dispersion Explorer</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
  body { font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 0; padding: 20px; background: #fafafa; color: #222; }
  h1 { font-size: 20px; margin-bottom: 4px; }
  .subtitle { color: #666; font-size: 13px; margin-bottom: 20px; }
  .controls { display: flex; gap: 16px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 16px; padding: 14px; background: white; border: 1px solid #ddd; border-radius: 6px; }
  .control { display: flex; flex-direction: column; gap: 4px; }
  .control label { font-size: 11px; text-transform: uppercase; color: #888; font-weight: 600; }
  select { padding: 6px 8px; border-radius: 4px; border: 1px solid #ccc; font-size: 13px; min-width: 160px; }
  .tabs { display: flex; gap: 4px; margin-bottom: 12px; }
  .tab { padding: 8px 16px; background: #eee; border-radius: 6px 6px 0 0; cursor: pointer; font-size: 13px; font-weight: 600; color: #555; }
  .tab.active { background: white; color: #111; border: 1px solid #ddd; border-bottom: none; }
  .panel { display: none; background: white; border: 1px solid #ddd; border-radius: 0 6px 6px 6px; padding: 16px; min-height: 500px; }
  .panel.active { display: block; }
  .note { font-size: 12px; color: #888; margin-top: 8px; }
  .msg { padding: 60px 20px; text-align: center; color: #888; font-size: 14px; }
  #plot { width: 100%; height: 640px; }
</style>
</head>
<body>

<h1>DGLM Dispersion Explorer</h1>
<div class="subtitle" id="subtitle"></div>

<div class="controls">
  <div class="control">
    <label>Model</label>
    <select id="modelSel"></select>
  </div>
  <div class="control">
    <label id="unitLabel">Cell type</label>
    <select id="unitSel"></select>
  </div>
  <div class="control">
    <label>Region</label>
    <select id="regionSel"></select>
  </div>
  <div class="control">
    <label>Significance (LFSR &lt;)</label>
    <select id="threshSel"></select>
  </div>
</div>

<div class="tabs">
  <div class="tab active" data-tab="methods">Methods</div>
  <div class="tab" data-tab="bar">Bar chart</div>
  <div class="tab" data-tab="volcano">Volcano</div>
  <div class="tab" data-tab="heatmap">Heatmap</div>
  <div class="tab" data-tab="gene">Gene search</div>
</div>

<div class="panel active" id="panel-methods" style="line-height:1.6;">
  <h2 style="margin-top:0;">What this is</h2>
  <p>DGLM dispersion analysis (mean ~ age+sex+mean_n_umi+n_cells, dispersion ~ age) on the macaque brain aging dataset, at louvain-subcluster x region resolution, with mashr run <b>separately per cell type</b> (each cell type's subclusters x regions form one mash fit) rather than pooled across all cell types at once -- a deliberate change from the prior single-pooled-mashr approach.</p>

  <h2>Pipeline rebuilt this round</h2>
  <p><b>Pseudobulk aggregation fixed:</b> previously averaged raw counts per animal instead of summing them, which broke CPM-based library-size normalization downstream. Now sums raw counts per animal x louvain-subcluster x region, with no per-animal cell-count minimum and no crude gene-presence filter at this stage.</p>
  <p><b>Gene filtering moved downstream and made explicit:</b> a gene is kept if its CPM is at least a chosen cutoff (0.5 used here) in at least 50% of that combo's animals -- computed from raw (non-TMM-normalized) library sizes, inclusive (&gt;=) at both the CPM and animal-count thresholds.</p>
  <p><b>Shat mode: sqrt(bvar) only.</b> Matches Chiou et al.'s original method exactly. The raw-bvar alternative was tested and found poorly calibrated (see prior methods notes) -- this build only offers the sqrt version, since that's the one actually used going forward.</p>

  <p><b>The Model dropdown above currently has a single option: sqrt-shat, louvain resolution, mashr run per cell type.</b> Earlier builds compared four combinations (raw/sqrt x celltype/louvain); this pipeline no longer produces the other three, so they're not offered.</p>
</div>

<div class="panel" id="panel-bar">
  <div class="controls" style="margin-bottom:12px;">
    <div class="control">
      <label>Data source</label>
      <select id="barSourceSel">
        <option value="mash">Mashr (mash_beta / mash_lfsr)</option>
        <option value="raw">Raw DGLM (beta / raw q&lt;0.05)</option>
      </select>
    </div>
  </div>
  <div id="plot-bar" style="width:100%;height:640px;"></div>
  <div class="note" id="note-bar"></div>
</div>

<div class="panel" id="panel-volcano">
  <div id="plot-volcano" style="width:100%;height:640px;"></div>
  <div class="note">mash_beta vs. mash_lfsr. Non-significant points are a random subsample (up to 2000) for file-size reasons; the full significant set is always shown, colored by direction. Hover a point to see its gene symbol and region. Pick a specific cell type (not "All") to view -- a volcano needs one unit at a time.</div>
</div>

<div class="panel" id="panel-heatmap">
  <div class="controls" style="margin-bottom:12px;">
    <div class="control">
      <label>Value</label>
      <select id="heatmapMetricSel"></select>
    </div>
  </div>
  <div id="plot-heatmap" style="width:100%;height:640px;"></div>
  <div class="note" id="note-heatmap"></div>
</div>

<div class="panel" id="panel-gene">
  <div class="controls" style="margin-bottom:12px;">
    <div class="control" style="min-width:280px;">
      <label>Gene symbol</label>
      <input id="geneSearchInput" type="text" placeholder="e.g. ADRB2" style="padding:6px 8px; border-radius:4px; border:1px solid #ccc; font-size:14px;">
    </div>
  </div>
  <div id="gene-result"></div>
  <div class="note">This searches the data already loaded for the currently selected Model and Cell type -- it shows mash_beta/mash_lfsr per region for that gene, not literal per-animal expression values (those aren't embedded here). Only significant genes (LFSR &lt; 0.5) and a random background sample are embedded per cell type, so a non-significant gene outside that sample won't be found here -- that doesn't mean it doesn't exist in the full data. For the real per-animal panel-b plot, use the command shown below once a gene is found.</div>
</div>

<script>
const THRESH_OPTIONS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5];
const DEFAULT_THRESH = 0.2;

const ABBR = {
  "astrocytes": "AST", "basket cells": "BC", "cerebellar neurons": "CER",
  "ependymal cells": "EPEN", "gabaergic neurons": "INH", "glutamatergic neurons": "EXC",
  "gabaergic": "INH", "glutamatergic": "EXC", "opc": "OPC", "opcs": "OPC",
  "medium spiny neurons": "MSN", "microglia": "MGL", "midbrain neurons": "MBN",
  "oligodendrocytes": "OLIG", "vascular cells": "VASC"
};
function normKey(c) {
  return String(c).trim().toLowerCase().replace(/[-_]/g, ' ');
}
function abbr(c) {
  const key = normKey(c);
  if (ABBR[key]) return ABBR[key];
  const m = String(c).match(/^([a-zA-Z_]+)_(\d+)$/);
  if (m) {
    const parentAbbr = ABBR[normKey(m[1])];
    if (parentAbbr) return parentAbbr + '-' + m[2];
  }
  return String(c);
}
function isort(arr) {
  return [...arr].sort((a, b) => String(a).toLowerCase().localeCompare(String(b).toLowerCase()));
}
function directionLabels(modelKey) {
  return { pos: 'Increase', neg: 'Decrease', posColor: '#d9534f', negColor: '#428bca' };
}
function metricLabel(key, modelKey) {
  const dl = directionLabels(modelKey);
  const overrides = {
    n_increase: '# sig genes, ' + dl.pos.toLowerCase(),
    n_decrease: '# sig genes, ' + dl.neg.toLowerCase(),
    n_net: 'Net (' + dl.pos + ' \u2212 ' + dl.neg + ')'
  };
  return overrides[key] || null;
}
function populateHeatmapMetric() {
  const m = currentModel();
  const sel = document.getElementById('heatmapMetricSel');
  const prev = sel.value;
  sel.innerHTML = '';
  const thresh = document.getElementById('threshSel').value;
  const byThresh = m.heatmap.by_threshold[thresh] || {};
  const allMetrics = { ...m.heatmap.static, ...byThresh };
  const preferredOrder = ['pct_sig', 'mash_beta_sig', 'mash_beta_mag_sig', 'n_increase', 'n_decrease',
                           'n_net', 'strongest_hit', 'genes_tested', 'raw_beta_all', 'raw_beta_sig'];
  preferredOrder.filter(k => allMetrics[k]).forEach(k => {
    const opt = document.createElement('option');
    opt.value = k;
    opt.textContent = metricLabel(k, modelSel.value) || allMetrics[k].label;
    sel.appendChild(opt);
  });
  if ([...sel.options].some(o => o.value === prev)) sel.value = prev;
}
function populateThresh() {
  const sel = document.getElementById('threshSel');
  const prev = sel.value;
  sel.innerHTML = '';
  THRESH_OPTIONS.forEach(t => {
    const opt = document.createElement('option');
    opt.value = String(t);
    opt.textContent = String(t);
    sel.appendChild(opt);
  });
  sel.value = [...sel.options].some(o => o.value === prev) ? prev : String(DEFAULT_THRESH);
}
function updateSubtitle() {
  document.getElementById('subtitle').textContent =
    'Transcriptional dispersion aging effects \u2014 mash-shrunk results (significance: mash LFSR < ' +
    document.getElementById('threshSel').value + ')';
}

const modelSel = document.getElementById('modelSel');
const unitSel = document.getElementById('unitSel');
const regionSel = document.getElementById('regionSel');
const unitLabel = document.getElementById('unitLabel');

Object.keys(DATA).forEach(k => {
  const opt = document.createElement('option');
  opt.value = k;
  opt.textContent = DATA[k].label;
  modelSel.appendChild(opt);
});

function currentModel() { return DATA[modelSel.value]; }

function populateUnitAndRegion() {
  const m = currentModel();
  unitLabel.textContent = (m.resolution === 'louvain') ? 'Subcluster' : 'Cell type';
  unitSel.innerHTML = '';
  const allOpt = document.createElement('option');
  allOpt.value = '__ALL__';
  allOpt.textContent = 'All (rollup)';
  unitSel.appendChild(allOpt);
  isort(m.units).forEach(u => {
    const opt = document.createElement('option');
    opt.value = u; opt.textContent = abbr(u) + (abbr(u) !== u ? ' (' + u + ')' : '');
    unitSel.appendChild(opt);
  });
  regionSel.innerHTML = '';
  const allRegionOpt = document.createElement('option');
  allRegionOpt.value = '__ALL__';
  allRegionOpt.textContent = 'All regions';
  regionSel.appendChild(allRegionOpt);
  isort(m.regions).forEach(r => {
    const opt = document.createElement('option');
    opt.value = r; opt.textContent = r;
    regionSel.appendChild(opt);
  });
}

function drawBar() {
  const m = currentModel();
  const unit = unitSel.value;
  const source = document.getElementById('barSourceSel').value;
  const dl = source === 'raw' ? { pos: 'Increase', neg: 'Decrease', posColor: '#d9534f', negColor: '#428bca' }
                               : directionLabels(modelSel.value);
  const thresh = document.getElementById('threshSel').value;
  const barData = source === 'raw' ? m.bar.raw : m.bar.mash[thresh];
  let regions, pos, neg;
  if (unit === '__ALL__') {
    regions = isort(Object.keys(barData.by_region));
    pos = regions.map(r => barData.by_region[r].increase);
    neg = regions.map(r => -barData.by_region[r].decrease);
  } else {
    const perRegion = barData.by_celltype_region[unit] || {};
    regions = isort(Object.keys(perRegion));
    pos = regions.map(r => perRegion[r].increase);
    neg = regions.map(r => -perRegion[r].decrease);
  }
  document.getElementById('note-bar').textContent = source === 'raw'
    ? 'Counts of genes significant at raw q<0.05, direction from raw beta sign \u2014 independent of the LFSR threshold selector above. Raw q-values are known to be poorly calibrated \u2014 interpret cautiously relative to the mashr view.'
    : ('Counts of significant genes (LFSR < ' + thresh + ') per region for the selected unit ("All" = rollup across all units). ' +
       'Positive = increases with age, negative = decreases.');
  Plotly.react('plot-bar', [
    { x: regions, y: pos, name: dl.pos, type: 'bar', marker: { color: dl.posColor } },
    { x: regions, y: neg, name: dl.neg, type: 'bar', marker: { color: dl.negColor } }
  ], {
    barmode: 'relative', title: (unit === '__ALL__' ? 'All units' : abbr(unit)) + ' \u2014 sig genes by region',
    yaxis: { title: 'Number of genes' }, margin: { t: 40 }
  }, {responsive: true});
}

function drawVolcano() {
  const m = currentModel();
  const dl = directionLabels(modelSel.value);
  const thresh = parseFloat(document.getElementById('threshSel').value);
  const container = document.getElementById('plot-volcano');
  if (unitSel.value === '__ALL__') {
    container.innerHTML = '<div class="msg">Pick a specific cell type (not "All") to view a volcano plot.</div>';
    return;
  }
  const unit = unitSel.value;
  let pts;
  if (regionSel.value === '__ALL__') {
    pts = [];
    Object.keys(m.volcano[unit] || {}).forEach(r => { pts = pts.concat(m.volcano[unit][r]); });
  } else {
    pts = (m.volcano[unit] && m.volcano[unit][regionSel.value]) || [];
  }
  const bgPts = pts.filter(p => p.l >= thresh);
  const posPts = pts.filter(p => p.l < thresh && p.b > 0);
  const negPts = pts.filter(p => p.l < thresh && p.b <= 0);
  const yOf = p => -Math.log10(Math.max(p.l, 1e-300));
  Plotly.newPlot('plot-volcano', [
    { x: bgPts.map(p => p.b), y: bgPts.map(yOf),
      text: bgPts.map(p => p.g), mode: 'markers', type: 'scattergl', name: 'Not sig (sample)',
      marker: { color: '#ccc', size: 5 }, hoverinfo: 'text+x+y' },
    { x: negPts.map(p => p.b), y: negPts.map(yOf),
      text: negPts.map(p => p.g), mode: 'markers', type: 'scattergl',
      name: dl.neg, marker: { color: dl.negColor, size: 7 }, hoverinfo: 'text+x+y' },
    { x: posPts.map(p => p.b), y: posPts.map(yOf),
      text: posPts.map(p => p.g), mode: 'markers', type: 'scattergl',
      name: dl.pos, marker: { color: dl.posColor, size: 7 }, hoverinfo: 'text+x+y' }
  ], {
    title: abbr(unit) + ' \u00d7 ' + (regionSel.value === '__ALL__' ? 'All regions' : regionSel.value),
    xaxis: { title: 'mash_beta' },
    yaxis: { title: '-log10(mash_lfsr)' }, margin: { t: 40 }
  }, {responsive: true});
}

function drawHeatmap() {
  const m = currentModel();
  const metricKey = document.getElementById('heatmapMetricSel').value;
  const thresh = document.getElementById('threshSel').value;
  const metric = m.heatmap.static[metricKey] || (m.heatmap.by_threshold[thresh] || {})[metricKey];
  if (!metric) return;

  const z = metric.z;
  const label = metricLabel(metricKey, modelSel.value) || metric.label;

  const flat = z.flat().filter(v => v !== null && v !== undefined);

  const yLabels = m.heatmap.units.map(abbr);

  let trace;
  if (metric.kind === 'sequential') {
    const maxV = flat.length ? Math.max(...flat) : 1;
    const seqColors = (metricKey === 'n_decrease')
      ? [[0, '#ffffff'], [1, '#2166ac']]
      : [[0, '#ffffff'], [1, '#b2182b']];
    trace = {
      z: z, x: m.heatmap.regions, y: m.heatmap.units, type: 'heatmap',
      colorscale: seqColors, zmin: 0, zmax: maxV || 1,
      colorbar: { title: label }
    };
  } else {
    const bound = flat.length ? Math.max(...flat.map(Math.abs)) : 1;
    trace = {
      z: z, x: m.heatmap.regions, y: m.heatmap.units, type: 'heatmap',
      colorscale: [[0, '#2166ac'], [0.5, '#ffffff'], [1, '#b2182b']],
      zmin: -bound, zmax: bound,
      colorbar: { title: label }
    };
  }
  trace.hovertemplate = '%{y} \u00d7 %{x}<br>' + label + ': %{z}<extra></extra>';

  Plotly.react('plot-heatmap', [trace], {
    title: label,
    margin: { t: 40, l: 140 },
    yaxis: { tickvals: m.heatmap.units, ticktext: yLabels, autorange: 'reversed' }
  }, {responsive: true});

  document.getElementById('note-heatmap').textContent =
    metricKey.startsWith('raw_beta')
      ? 'Raw per-gene q-values are known to be poorly calibrated \u2014 interpret this metric cautiously relative to the mash-based ones. Region selector doesn\'t apply on this tab; this metric doesn\'t depend on the LFSR threshold.'
      : 'Region selector doesn\'t apply on this tab \u2014 heatmap always shows all regions. Blank cells = that unit x region combination doesn\'t exist in the data. Uses the selected LFSR threshold (' + thresh + ').';
}

function showUnavailable() {
  const msg = '<div class="msg">This model/resolution hasn\'t been (re)generated yet \u2014 check back once the pipeline finishes.</div>';
  document.getElementById('plot-bar').innerHTML = msg;
  document.getElementById('note-bar').textContent = '';
  document.getElementById('plot-volcano').innerHTML = msg;
  document.getElementById('plot-heatmap').innerHTML = msg;
  document.getElementById('note-heatmap').textContent = '';
  document.getElementById('gene-result').innerHTML = msg;
}

const CHECKPOINTS_DIR = {
  v3_sqrt_louvain: '/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'
};

function searchGene() {
  const container = document.getElementById('gene-result');
  const query = document.getElementById('geneSearchInput').value.trim();
  if (!query) {
    container.innerHTML = '<div class="msg">Type a gene symbol above to search.</div>';
    return;
  }
  const m = currentModel();
  if (!m.available) { showUnavailable(); return; }
  const unit = unitSel.value;
  if (unit === '__ALL__') {
    container.innerHTML = '<div class="msg">Pick a specific cell type (not "All") to search a gene within it.</div>';
    return;
  }
  const byRegion = m.volcano[unit] || {};
  const hits = [];
  Object.keys(byRegion).forEach(region => {
    byRegion[region].forEach(p => {
      if (p.g.toLowerCase() === query.toLowerCase()) hits.push({ region, b: p.b, l: p.l });
    });
  });

  container.innerHTML = '';
  if (hits.length === 0) {
    const msgDiv = document.createElement('div');
    msgDiv.className = 'msg';
    msgDiv.textContent = '"' + query + '" not found in the visualized sample for ' + abbr(unit) +
      ' in this model. Either it\'s not significant here, or it fell outside the random background sample embedded for file-size reasons. It may still exist in the full data -- the command below will check.';
    container.appendChild(msgDiv);
  } else {
    hits.sort((a, b) => a.region.localeCompare(b.region));
    const div = document.createElement('div');
    div.id = 'plot-gene';
    div.style.height = '400px';
    container.appendChild(div);
    Plotly.newPlot('plot-gene', [{
      x: hits.map(h => h.region), y: hits.map(h => h.b), type: 'bar',
      marker: { color: hits.map(h => h.b > 0 ? '#d9534f' : '#428bca') },
      text: hits.map(h => 'lfsr=' + h.l.toExponential(2)), hoverinfo: 'x+y+text'
    }], {
      title: query + ' \u2014 ' + abbr(unit) + ', mash_beta by region',
      yaxis: { title: 'mash_beta' }, margin: { t: 40 }
    }, { responsive: true });
  }

  const ctDir = CHECKPOINTS_DIR[modelSel.value] + '/' + unit.replace(/_\d+$/, '') + '/dglm_checkpoints_cutoff0.5';
  const cmdBlock = document.createElement('div');
  cmdBlock.style.marginTop = '16px';
  cmdBlock.innerHTML = '<div class="note" style="margin-bottom:4px;">For the real per-animal fig4b plot for this exact gene and cell type:</div>' +
    '<pre style="background:#f4f4f4; padding:10px; border-radius:4px; font-size:12px; overflow-x:auto;">' +
    '/home/easmit31/.conda/envs/mashr_env/bin/Rscript /scratch/easmit31/dispersion/dglm/scripts/dglm_plot_fig4.R \\\n' +
    '  --checkpoints ' + ctDir + ' \\\n' +
    '  --master_tsv ' + ctDir + '/master_dglm_combined.tsv \\\n' +
    '  --figdir ' + ctDir + '/figures_gene_lookup \\\n' +
    '  --panel_b_only --example_symbol ' + query + ' --example_ct ' + unit +
    '</pre>';
  container.appendChild(cmdBlock);
}

function drawActive() {
  const m = currentModel();
  if (!m.available) { showUnavailable(); return; }
  const active = document.querySelector('.tab.active').dataset.tab;
  if (active === 'bar') drawBar();
  else if (active === 'volcano') drawVolcano();
  else if (active === 'heatmap') drawHeatmap();
  else if (active === 'gene') searchGene();
}

document.querySelectorAll('.tab').forEach(tab => {
  tab.addEventListener('click', () => {
    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
    tab.classList.add('active');
    document.getElementById('panel-' + tab.dataset.tab).classList.add('active');
    drawActive();
  });
});

modelSel.addEventListener('change', () => { populateUnitAndRegion(); populateHeatmapMetric(); drawActive(); });
unitSel.addEventListener('change', drawActive);
regionSel.addEventListener('change', drawActive);
document.getElementById('heatmapMetricSel').addEventListener('change', drawActive);
document.getElementById('barSourceSel').addEventListener('change', drawActive);
document.getElementById('threshSel').addEventListener('change', () => { updateSubtitle(); populateHeatmapMetric(); drawActive(); });
document.getElementById('geneSearchInput').addEventListener('input', () => {
  if (document.querySelector('.tab.active').dataset.tab === 'gene') searchGene();
});

populateThresh();
updateSubtitle();
populateUnitAndRegion();
populateHeatmapMetric();
drawActive();
</script>
</body>
</html>
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True, help='JSON file from build_dglm_explorer_data.R')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    with open(args.data) as f:
        data_json_text = f.read()

    # Split template at the DATA-consuming line, inject the real DATA object
    marker = "const modelSel = document.getElementById('modelSel');"
    idx = TEMPLATE_HEAD.find(marker)
    if idx == -1:
        raise SystemExit("Could not find injection point in template")

    html = TEMPLATE_HEAD[:idx] + f"const DATA = {data_json_text};\n\n" + TEMPLATE_HEAD[idx:]

    with open(args.out, 'w') as f:
        f.write(html)
    print(f"Saved: {args.out}")

if __name__ == '__main__':
    main()
