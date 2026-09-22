#!/usr/bin/env python3
"""
Single-file DGLM explorer builder. Reads master_dglm_*.tsv files directly via
pandas -- no R, no mashr, no .rds. Unit selector has three tiers: All, one
parent cell type (grouped/summed on bar; grouped/mean/median on heatmap;
grouped on volcano), or one specific subcluster. Natural sort throughout.
Heatmap: value-label toggle, mean/median aggregation, no stray gridlines.
A region checkbox filter (top controls bar) restricts which regions appear
across bar/heatmap/volcano, with plots autoscaling to whatever remains.

Usage:
    python build_dglm_explorer.py --out /scratch/easmit31/dispersion/dglm/dglm_explorer.html
"""
import argparse
import glob
import json
import math
import os
import re

import pandas as pd

BASE_DIR = '/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'
CUTOFF = '0.5'
THRESH_OPTIONS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
MAX_BACKGROUND = 2000

def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]

def nsorted(iterable):
    return sorted(iterable, key=natural_key)

def nz(v):
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    return v

def load_master(ckpt_suffix, tsv_name):
    frames = []
    cell_type_dirs = sorted(d for d in glob.glob(os.path.join(BASE_DIR, '*')) if os.path.isdir(d))
    for d in cell_type_dirs:
        tsv = os.path.join(d, f'dglm_checkpoints_cutoff{CUTOFF}{ckpt_suffix}', tsv_name)
        if not os.path.exists(tsv):
            continue
        df = pd.read_csv(tsv, sep='\t')
        parent = os.path.basename(d)
        df['parent_cell_type'] = parent
        is_opc_olig = df['cell_type'].astype(str).str.startswith('opc-olig_')
        if is_opc_olig.any():
            suffix = df.loc[is_opc_olig, 'cell_type'].str.replace('^opc-olig_', '', regex=True)
            df.loc[is_opc_olig, 'cell_type'] = parent + '_' + suffix
        frames.append(df)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)

def build_bar(df):
    mash_by_thresh = {}
    for thresh in THRESH_OPTIONS:
        sig = df['mash_lfsr'] < thresh
        grp = df[sig].groupby(['cell_type', 'region'])['mash_beta']
        inc = grp.apply(lambda s: int((s > 0).sum()))
        dec = grp.apply(lambda s: int((s <= 0).sum()))
        by_ct_region = {}
        for unit in nsorted(df['cell_type'].unique()):
            by_ct_region[unit] = {}
            for region in nsorted(df['region'].unique()):
                i = int(inc.get((unit, region), 0))
                de = int(dec.get((unit, region), 0))
                by_ct_region[unit][region] = {'increase': i, 'decrease': de}
        by_region = {}
        for region in nsorted(df['region'].unique()):
            i = sum(by_ct_region[u].get(region, {}).get('increase', 0) for u in by_ct_region)
            de = sum(by_ct_region[u].get(region, {}).get('decrease', 0) for u in by_ct_region)
            by_region[region] = {'increase': i, 'decrease': de}
        mash_by_thresh[str(thresh)] = {'by_region': by_region, 'by_celltype_region': by_ct_region}

    rsig = df['qvalue'] < 0.05
    rgrp = df[rsig].groupby(['cell_type', 'region'])['beta']
    rinc = rgrp.apply(lambda s: int((s > 0).sum()))
    rdec = rgrp.apply(lambda s: int((s <= 0).sum()))
    raw_by_ct_region = {}
    for unit in nsorted(df['cell_type'].unique()):
        raw_by_ct_region[unit] = {}
        for region in nsorted(df['region'].unique()):
            raw_by_ct_region[unit][region] = {
                'increase': int(rinc.get((unit, region), 0)),
                'decrease': int(rdec.get((unit, region), 0))
            }
    raw_by_region = {}
    for region in nsorted(df['region'].unique()):
        i = sum(raw_by_ct_region[u].get(region, {}).get('increase', 0) for u in raw_by_ct_region)
        de = sum(raw_by_ct_region[u].get(region, {}).get('decrease', 0) for u in raw_by_ct_region)
        raw_by_region[region] = {'increase': i, 'decrease': de}

    return {'mash': mash_by_thresh, 'raw': {'by_region': raw_by_region, 'by_celltype_region': raw_by_ct_region}}

def build_volcano(df):
    volcano = {}
    for unit, gdf in df.groupby('cell_type'):
        volcano[unit] = {}
        for region, rdf in gdf.groupby('region'):
            pts = []
            sig_rows = rdf[rdf['mash_lfsr'] < 0.5]
            bg_rows = rdf[rdf['mash_lfsr'] >= 0.5]
            if len(bg_rows) > MAX_BACKGROUND:
                bg_rows = bg_rows.sample(n=MAX_BACKGROUND, random_state=0)
            for _, row in pd.concat([sig_rows, bg_rows]).iterrows():
                sym = row['symbol'] if isinstance(row['symbol'], str) and row['symbol'] else row['ensembl_id']
                pts.append({'g': sym, 'b': nz(float(row['mash_beta'])), 'l': nz(float(row['mash_lfsr']))})
            volcano[unit][region] = pts
    return volcano

def build_heatmap(df):
    all_units = nsorted(df['cell_type'].unique())
    all_regions = nsorted(df['region'].unique())

    def mk_matrix(fill_fn):
        return [[fill_fn(u, r) for r in all_regions] for u in all_units]

    grouped = {k: v for k, v in df.groupby(['cell_type', 'region'])}

    def sub(u, r):
        return grouped.get((u, r))

    genes_tested = mk_matrix(lambda u, r: (len(sub(u, r)) if sub(u, r) is not None else None))
    raw_beta_all = mk_matrix(lambda u, r: (nz(float(sub(u, r)['beta'].abs().mean())) if sub(u, r) is not None else None))
    def raw_beta_sig_fn(u, r):
        s = sub(u, r)
        if s is None:
            return None
        sig = s[s['qvalue'] < 0.05]
        return nz(float(sig['beta'].abs().mean())) if len(sig) else 0
    raw_beta_sig = mk_matrix(raw_beta_sig_fn)

    by_threshold = {}
    for thresh in THRESH_OPTIONS:
        def pct_sig_fn(u, r, thresh=thresh):
            s = sub(u, r)
            if s is None:
                return None
            return 100 * float((s['mash_lfsr'] < thresh).sum()) / len(s)
        def n_inc_fn(u, r, thresh=thresh):
            s = sub(u, r)
            if s is None:
                return None
            return int(((s['mash_lfsr'] < thresh) & (s['mash_beta'] > 0)).sum())
        def n_dec_fn(u, r, thresh=thresh):
            s = sub(u, r)
            if s is None:
                return None
            return int(((s['mash_lfsr'] < thresh) & (s['mash_beta'] <= 0)).sum())
        def mash_beta_sig_fn(u, r, thresh=thresh):
            s = sub(u, r)
            if s is None:
                return None
            sig = s[s['mash_lfsr'] < thresh]
            return nz(float(sig['mash_beta'].mean())) if len(sig) else 0
        def mash_beta_mag_sig_fn(u, r, thresh=thresh):
            s = sub(u, r)
            if s is None:
                return None
            sig = s[s['mash_lfsr'] < thresh]
            return nz(float(sig['mash_beta'].abs().mean())) if len(sig) else 0

        n_increase = mk_matrix(n_inc_fn)
        n_decrease = mk_matrix(n_dec_fn)
        n_net = [[(a - b) if (a is not None and b is not None) else None
                  for a, b in zip(row_i, row_d)] for row_i, row_d in zip(n_increase, n_decrease)]

        by_threshold[str(thresh)] = {
            'pct_sig': {'z': mk_matrix(pct_sig_fn), 'kind': 'sequential', 'label': '% significant'},
            'n_increase': {'z': n_increase, 'kind': 'sequential', 'label': '# sig genes, increase'},
            'n_decrease': {'z': n_decrease, 'kind': 'sequential', 'label': '# sig genes, decrease'},
            'n_net': {'z': n_net, 'kind': 'diverging', 'label': 'Net (increase - decrease)'},
            'mash_beta_sig': {'z': mk_matrix(mash_beta_sig_fn), 'kind': 'diverging', 'label': 'Mean mash beta (sig genes, averaged)'},
            'mash_beta_mag_sig': {'z': mk_matrix(mash_beta_mag_sig_fn), 'kind': 'sequential', 'label': 'Mean |mash beta| (sig genes, averaged)'},
        }

    return {
        'units': all_units, 'regions': all_regions,
        'static': {
            'genes_tested': {'z': genes_tested, 'kind': 'sequential', 'label': 'Genes tested'},
            'raw_beta_all': {'z': raw_beta_all, 'kind': 'sequential', 'label': 'Mean |raw beta| (all genes, averaged)'},
            'raw_beta_sig': {'z': raw_beta_sig, 'kind': 'sequential', 'label': 'Mean |raw beta| (raw q<0.05, averaged)'},
        },
        'by_threshold': by_threshold,
    }

def build_model(ckpt_suffix, tsv_name, model_label):
    df = load_master(ckpt_suffix, tsv_name)
    if df is None:
        print(f'[{model_label}] no master TSVs found')
        return {'label': model_label, 'resolution': 'louvain', 'available': False,
                'units': [], 'regions': [], 'parents': [], 'unit_to_parent': {},
                'bar': {}, 'volcano': {}, 'heatmap': {}}
    all_units = nsorted(df['cell_type'].unique())
    all_regions = nsorted(df['region'].unique())
    unit_to_parent = df.drop_duplicates('cell_type').set_index('cell_type')['parent_cell_type'].to_dict()
    parents = nsorted(df['parent_cell_type'].unique())
    print(f'[{model_label}] cell types (subclusters): {len(all_units)} | parents: {len(parents)} | regions: {", ".join(all_regions)}')
    return {
        'label': model_label, 'resolution': 'louvain', 'available': True,
        'units': all_units, 'regions': all_regions,
        'parents': parents, 'unit_to_parent': unit_to_parent,
        'bar': build_bar(df), 'volcano': build_volcano(df), 'heatmap': build_heatmap(df),
    }

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
  .checkctrl { display: flex; align-items: center; gap: 6px; }
  .checkctrl input { width: 16px; height: 16px; }
  .regionBtn { padding: 6px 8px; border-radius: 4px; border: 1px solid #ccc; font-size: 13px; min-width: 160px; background: white; cursor: pointer; text-align: left; }
  .regionPanel { display: none; flex-wrap: wrap; gap: 6px 16px; padding: 12px 14px; background: white; border: 1px solid #ddd; border-radius: 6px; margin: -10px 0 16px 0; }
  .regionPanel label { display: flex; align-items: center; gap: 5px; font-size: 13px; cursor: pointer; text-transform: none; color: #333; }
  .regionPanel input { width: 15px; height: 15px; }
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
    <label>Region (single, volcano)</label>
    <select id="regionSel"></select>
  </div>
  <div class="control">
    <label>Regions shown (all tabs)</label>
    <button type="button" class="regionBtn" id="regionFilterToggle">All regions \u25be</button>
  </div>
  <div class="control">
    <label>Significance (LFSR &lt;)</label>
    <select id="threshSel"></select>
  </div>
</div>

<div class="regionPanel" id="regionFilterPanel"></div>

<div class="tabs">
  <div class="tab active" data-tab="methods">Methods</div>
  <div class="tab" data-tab="bar">Bar chart</div>
  <div class="tab" data-tab="volcano">Volcano</div>
  <div class="tab" data-tab="heatmap">Heatmap</div>
  <div class="tab" data-tab="gene">Gene search</div>
</div>

<div class="panel active" id="panel-methods" style="line-height:1.6;">
  <h2 style="margin-top:0;">What this shows</h2>
  <p>Age-associated changes in per-gene transcriptional variability (dispersion) within each cell type and brain region, in the macaque brain aging dataset. For each gene, this looks at whether the <i>spread</i> of its expression across single cells changes with age -- not just the average level. A gene can show no change in mean expression but still become more or less variable with age, which is what this analysis is built to detect.</p>
  <p>Two models are available from the Model selector above: dispersion as a function of age alone, and dispersion as a function of age and sex together (sex included as a covariate).</p>

  <h2>The data</h2>
  <p>Single-nucleus RNA-seq from multiple brain regions of rhesus macaques spanning a range of adult ages. Cells within each broad cell type (astrocytes, microglia, glutamatergic neurons, etc.) are further split into louvain subclusters -- finer-grained groupings of transcriptionally similar cells within that cell type. Every plot here can be viewed at either the subcluster level or rolled up to the whole cell type.</p>

  <h2>The model</h2>
  <p>Each gene is fit with a double generalized linear model (DGLM): one submodel for the mean level of expression (mean ~ age + sex + sequencing depth + cells per sample), and a separate submodel for the dispersion, or spread, of expression around that mean (dispersion ~ age, or dispersion ~ age + sex). Fitting these separately is what lets the analysis distinguish "this gene is more/less variable with age" from "this gene's average level changed with age" -- the two aren't the same thing, and a gene can show one without the other.</p>
  <p>Before a gene is tested in a given subcluster x region, it has to clear a minimum expression filter there, applied independently per subcluster x region. Because of that, the number of genes tested (and so how many can turn up significant) varies across cell types and regions -- rarer or smaller subclusters test fewer genes than large, abundant ones. That's expected, not a data quality issue.</p>

  <h2>Statistical pooling</h2>
  <p>Raw per-gene effect estimates are noisy, especially for smaller subclusters. To get more reliable effect sizes and significance calls, results from all of a cell type's own subclusters and regions are jointly modeled together using mashr (multivariate adaptive shrinkage) -- pooling information across a cell type's own conditions to shrink noisy estimates toward more plausible values. This is done separately per cell type, not pooled across all cell types at once.</p>

  <h2>How to read the tabs</h2>
  <p><b>Bar chart</b> -- number of significant genes per region, split by whether dispersion increases or decreases with age.<br>
  <b>Volcano</b> -- effect size vs. significance for individual genes, for one cell type or subcluster at a time.<br>
  <b>Heatmap</b> -- a summary metric (gene counts, % significant, or average effect size) across every cell type/subcluster and region at once.<br>
  <b>Gene search</b> -- look up how a specific gene behaves across regions, for the selected cell type.</p>
  <p>The "Regions shown" control in the top bar restricts which regions appear across every tab -- uncheck a region and the plots redraw and autoscale to whatever remains checked.</p>

  <h2>On significance</h2>
  <p>Two significance measures appear throughout: mash LFSR (local false sign rate -- the primary one, threshold adjustable above) and raw DGLM q-values. LFSR is the probability that the estimated <i>direction</i> of an effect (increase vs. decrease) is wrong -- a different quantity than a p-value, and low LFSR reflects confidence in the direction, not necessarily a large effect size. The raw q-values are known to be poorly calibrated on this dataset and should be interpreted cautiously; mash LFSR is the more trustworthy of the two.</p>
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
    <div class="control" id="barGroupModeControl" style="display:none;">
      <label>Cell type view</label>
      <select id="barGroupModeSel">
        <option value="grouped">Grouped by subcluster (small charts)</option>
        <option value="summed">Summed across its subclusters</option>
      </select>
    </div>
  </div>
  <div id="plot-bar" style="width:100%;height:640px;"></div>
  <div class="note" id="note-bar"></div>
</div>

<div class="panel" id="panel-volcano">
  <div id="plot-volcano" style="width:100%;height:640px;"></div>
  <div class="note">mash_beta vs. mash_lfsr. Non-significant points are a random subsample (up to 2000) for file-size reasons; the full significant set is always shown, colored by direction. Hover a point to see its gene symbol and region.</div>
</div>

<div class="panel" id="panel-heatmap">
  <div class="controls" style="margin-bottom:12px;">
    <div class="control">
      <label>Value</label>
      <select id="heatmapMetricSel"></select>
    </div>
    <div class="control" id="heatmapAggControl">
      <label>Rows</label>
      <select id="heatmapAggSel">
        <option value="none">Every subcluster</option>
        <option value="mean">Mean per cell type</option>
        <option value="median">Median per cell type</option>
      </select>
    </div>
    <div class="control">
      <label>&nbsp;</label>
      <div class="checkctrl">
        <input type="checkbox" id="heatmapShowValues">
        <label for="heatmapShowValues" style="text-transform:none; font-size:13px; color:#333;">Show values on cells</label>
      </div>
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
  <div class="note">Searches the data loaded for the currently selected Model and unit (must be a single subcluster, not "All" or a grouped cell type). Only significant genes (LFSR &lt; 0.5) and a random background sample are embedded per subcluster, so a non-significant gene outside that sample won't be found here even though it may exist in the full data.</div>
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
  const chunk = s => (String(s).match(/(\d+|\D+)/g) || []);
  return [...arr].sort((a, b) => {
    const ac = chunk(a), bc = chunk(b);
    const len = Math.max(ac.length, bc.length);
    for (let i = 0; i < len; i++) {
      const ap = ac[i] || '', bp = bc[i] || '';
      const aNum = /^\d+$/.test(ap), bNum = /^\d+$/.test(bp);
      if (aNum && bNum) {
        const diff = parseInt(ap, 10) - parseInt(bp, 10);
        if (diff !== 0) return diff;
      } else {
        const cmp = ap.toLowerCase().localeCompare(bp.toLowerCase());
        if (cmp !== 0) return cmp;
      }
    }
    return 0;
  });
}
function directionLabels() {
  return { pos: 'Increase', neg: 'Decrease', posColor: '#d9534f', negColor: '#428bca' };
}
function metricLabel(key) {
  const dl = directionLabels();
  const overrides = {
    n_increase: '# sig genes, ' + dl.pos.toLowerCase(),
    n_decrease: '# sig genes, ' + dl.neg.toLowerCase(),
    n_net: 'Net (' + dl.pos + ' \u2212 ' + dl.neg + ')'
  };
  return overrides[key] || null;
}
function fmtVal(v) {
  if (v === null || v === undefined) return '';
  if (Number.isInteger(v)) return String(v);
  return v.toFixed(2);
}
function median(vals) {
  const s = [...vals].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 !== 0 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
}
function aggregateRows(units, z, groupOf, mode) {
  const groups = isort([...new Set(units.map(groupOf))]);
  const nCols = z[0] ? z[0].length : 0;
  const aggZ = groups.map(g => {
    const rowIdxs = units.map((u, i) => groupOf(u) === g ? i : -1).filter(i => i >= 0);
    const outRow = [];
    for (let c = 0; c < nCols; c++) {
      const vals = rowIdxs.map(i => z[i][c]).filter(v => v !== null && v !== undefined);
      if (vals.length === 0) { outRow.push(null); continue; }
      outRow.push(mode === 'median' ? median(vals) : vals.reduce((a, b) => a + b, 0) / vals.length);
    }
    return outRow;
  });
  return { units: groups, z: aggZ };
}
function checkedRegions() {
  return [...document.querySelectorAll('.regionFilterCb:checked')].map(cb => cb.value);
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
                           'n_net', 'genes_tested', 'raw_beta_all', 'raw_beta_sig'];
  preferredOrder.filter(k => allMetrics[k]).forEach(k => {
    const opt = document.createElement('option');
    opt.value = k;
    opt.textContent = metricLabel(k) || allMetrics[k].label;
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

function subclustersOfParent(m, parent) {
  return isort(m.units.filter(u => m.unit_to_parent[u] === parent));
}

function updateControlVisibility() {
  const unit = unitSel.value;
  document.getElementById('barGroupModeControl').style.display = unit.startsWith('PARENT:') ? 'flex' : 'none';
  document.getElementById('heatmapAggControl').style.display =
    (unit === '__ALL__' || unit.startsWith('PARENT:')) ? 'flex' : 'none';
}

function updateRegionFilterButtonLabel() {
  const all = document.querySelectorAll('.regionFilterCb').length;
  const checked = checkedRegions().length;
  document.getElementById('regionFilterToggle').textContent =
    (checked === all ? 'All regions' : checked + ' of ' + all + ' regions') + ' \u25be';
}

function populateRegionSelOptions() {
  const checked = isort(checkedRegions());
  const prev = regionSel.value;
  regionSel.innerHTML = '';
  const allRegionOpt = document.createElement('option');
  allRegionOpt.value = '__ALL__';
  allRegionOpt.textContent = 'All regions';
  regionSel.appendChild(allRegionOpt);
  checked.forEach(r => {
    const opt = document.createElement('option');
    opt.value = r; opt.textContent = r;
    regionSel.appendChild(opt);
  });
  regionSel.value = [...regionSel.options].some(o => o.value === prev) ? prev : '__ALL__';
}

function onRegionFilterChange() {
  updateRegionFilterButtonLabel();
  populateRegionSelOptions();
  drawActive();
}

function populateRegionFilter() {
  const m = currentModel();
  const panel = document.getElementById('regionFilterPanel');
  panel.innerHTML = '';
  isort(m.regions).forEach(r => {
    const label = document.createElement('label');
    const cb = document.createElement('input');
    cb.type = 'checkbox'; cb.value = r; cb.checked = true; cb.className = 'regionFilterCb';
    cb.addEventListener('change', onRegionFilterChange);
    label.appendChild(cb);
    label.appendChild(document.createTextNode(r));
    panel.appendChild(label);
  });
  updateRegionFilterButtonLabel();
}

document.getElementById('regionFilterToggle').addEventListener('click', () => {
  const panel = document.getElementById('regionFilterPanel');
  panel.style.display = panel.style.display === 'flex' ? 'none' : 'flex';
});

function populateUnitAndRegion() {
  const m = currentModel();
  unitLabel.textContent = 'Cell type / subcluster';
  unitSel.innerHTML = '';

  const allOpt = document.createElement('option');
  allOpt.value = '__ALL__';
  allOpt.textContent = 'All (summed across every subcluster)';
  unitSel.appendChild(allOpt);

  const parentGroup = document.createElement('optgroup');
  parentGroup.label = 'Cell type (bar: grouped or summed; heatmap: grouped, mean, or median; volcano: grouped)';
  isort(m.parents).forEach(p => {
    const opt = document.createElement('option');
    opt.value = 'PARENT:' + p;
    opt.textContent = abbr(p) + ' (' + p + ')';
    parentGroup.appendChild(opt);
  });
  unitSel.appendChild(parentGroup);

  const subGroup = document.createElement('optgroup');
  subGroup.label = 'Individual subclusters';
  isort(m.units).forEach(u => {
    const opt = document.createElement('option');
    opt.value = u; opt.textContent = abbr(u) + (abbr(u) !== u ? ' (' + u + ')' : '');
    subGroup.appendChild(opt);
  });
  unitSel.appendChild(subGroup);

  populateRegionFilter();
  populateRegionSelOptions();
  updateControlVisibility();
}

function gridDims(n) {
  const cols = Math.min(4, Math.max(1, Math.ceil(Math.sqrt(n))));
  const rows = Math.ceil(n / cols);
  return { rows, cols };
}

function sizeContainerForGrid(container, rows) {
  const h = Math.max(640, rows * 260);
  container.style.height = h + 'px';
  return h;
}

function drawBar() {
  const m = currentModel();
  const unit = unitSel.value;
  const source = document.getElementById('barSourceSel').value;
  const dl = directionLabels();
  const thresh = document.getElementById('threshSel').value;
  const barData = source === 'raw' ? m.bar.raw : m.bar.mash[thresh];
  const noteEl = document.getElementById('note-bar');
  const container = document.getElementById('plot-bar');
  const checked = checkedRegions();

  if (unit.startsWith('PARENT:')) {
    const parent = unit.slice('PARENT:'.length);
    const subs = subclustersOfParent(m, parent);
    const mode = document.getElementById('barGroupModeSel').value;

    if (mode === 'summed') {
      container.style.height = '640px';
      const regionSet = new Set();
      subs.forEach(s => Object.keys(barData.by_celltype_region[s] || {}).forEach(r => regionSet.add(r)));
      const regions = isort([...regionSet].filter(r => checked.includes(r)));
      const pos = regions.map(r => subs.reduce((acc, s) => acc + ((barData.by_celltype_region[s] || {})[r]?.increase || 0), 0));
      const neg = regions.map(r => -subs.reduce((acc, s) => acc + ((barData.by_celltype_region[s] || {})[r]?.decrease || 0), 0));
      noteEl.textContent = 'Summed across all ' + subs.length + ' subclusters of ' + abbr(parent) + ' (' + parent + ').';
      Plotly.newPlot(container, [
        { x: regions, y: pos, name: dl.pos, type: 'bar', marker: { color: dl.posColor } },
        { x: regions, y: neg, name: dl.neg, type: 'bar', marker: { color: dl.negColor } }
      ], {
        barmode: 'relative', title: abbr(parent) + ' (summed across subclusters) \u2014 sig genes by region',
        yaxis: { title: 'Number of genes' }, margin: { t: 40 }
      }, {responsive: true});
      return;
    }

    const { rows, cols } = gridDims(subs.length);
    sizeContainerForGrid(container, rows);
    const traces = [];
    const annotations = [];
    subs.forEach((sub, i) => {
      const axisNum = i === 0 ? '' : String(i + 1);
      const perRegion = barData.by_celltype_region[sub] || {};
      const regions = isort(Object.keys(perRegion).filter(r => checked.includes(r)));
      const pos = regions.map(r => perRegion[r].increase);
      const neg = regions.map(r => -perRegion[r].decrease);
      traces.push({ x: regions, y: pos, name: dl.pos, type: 'bar', marker: { color: dl.posColor },
                    xaxis: 'x' + axisNum, yaxis: 'y' + axisNum, showlegend: i === 0 });
      traces.push({ x: regions, y: neg, name: dl.neg, type: 'bar', marker: { color: dl.negColor },
                    xaxis: 'x' + axisNum, yaxis: 'y' + axisNum, showlegend: i === 0 });
      annotations.push({ text: abbr(sub), showarrow: false, xref: 'x' + axisNum + ' domain', yref: 'y' + axisNum + ' domain',
                          x: 0.5, y: 1.18, font: { size: 12 } });
    });
    noteEl.textContent = 'One small chart per subcluster of ' + abbr(parent) + ' (' + parent + ') \u2014 each shows that subcluster\'s own sig-gene counts by region, not summed together. Switch to "Summed" above to see them combined.';
    Plotly.newPlot(container, traces, {
      grid: { rows, columns: cols, pattern: 'independent', xgap: 0.35, ygap: 0.5 },
      barmode: 'relative', title: abbr(parent) + ' \u2014 all subclusters',
      annotations, margin: { t: 60, b: 40 }, showlegend: true,
      height: sizeContainerForGrid(container, rows)
    }, { responsive: true });
    return;
  }

  container.style.height = '640px';
  let regions, pos, neg;
  if (unit === '__ALL__') {
    regions = isort(Object.keys(barData.by_region).filter(r => checked.includes(r)));
    pos = regions.map(r => barData.by_region[r].increase);
    neg = regions.map(r => -barData.by_region[r].decrease);
  } else {
    const perRegion = barData.by_celltype_region[unit] || {};
    regions = isort(Object.keys(perRegion).filter(r => checked.includes(r)));
    pos = regions.map(r => perRegion[r].increase);
    neg = regions.map(r => -perRegion[r].decrease);
  }
  noteEl.textContent = source === 'raw'
    ? 'Counts of genes significant at raw q<0.05, direction from raw beta sign \u2014 independent of the LFSR threshold selector above. Raw q-values are known to be poorly calibrated \u2014 interpret cautiously relative to the mashr view.'
    : ('Counts of significant genes (LFSR < ' + thresh + ') per region for the selected unit (' +
       (unit === '__ALL__' ? 'summed across every subcluster' : 'this subcluster only') + '). ' +
       'Positive = increases with age, negative = decreases.');
  Plotly.newPlot(container, [
    { x: regions, y: pos, name: dl.pos, type: 'bar', marker: { color: dl.posColor } },
    { x: regions, y: neg, name: dl.neg, type: 'bar', marker: { color: dl.negColor } }
  ], {
    barmode: 'relative', title: (unit === '__ALL__' ? 'All units (summed)' : abbr(unit)) + ' \u2014 sig genes by region',
    yaxis: { title: 'Number of genes' }, margin: { t: 40 }
  }, {responsive: true});
}

function volcanoTracesFor(m, unit, dl, thresh) {
  let pts;
  if (regionSel.value === '__ALL__') {
    const checked = checkedRegions();
    pts = [];
    Object.keys(m.volcano[unit] || {}).forEach(r => { if (checked.includes(r)) pts = pts.concat(m.volcano[unit][r]); });
  } else {
    pts = (m.volcano[unit] && m.volcano[unit][regionSel.value]) || [];
  }
  const bgPts = pts.filter(p => p.l >= thresh);
  const posPts = pts.filter(p => p.l < thresh && p.b > 0);
  const negPts = pts.filter(p => p.l < thresh && p.b <= 0);
  const yOf = p => -Math.log10(Math.max(p.l, 1e-300));
  return { bgPts, posPts, negPts, yOf };
}

function drawVolcano() {
  const m = currentModel();
  const dl = directionLabels();
  const thresh = parseFloat(document.getElementById('threshSel').value);
  const container = document.getElementById('plot-volcano');

  if (unitSel.value === '__ALL__') {
    container.style.height = '640px';
    container.innerHTML = '<div class="msg">Pick a cell type or a specific subcluster (not "All") to view a volcano plot.</div>';
    return;
  }

  if (unitSel.value.startsWith('PARENT:')) {
    const parent = unitSel.value.slice('PARENT:'.length);
    const subs = subclustersOfParent(m, parent);
    const { rows, cols } = gridDims(subs.length);
    sizeContainerForGrid(container, rows);
    const traces = [];
    const annotations = [];
    subs.forEach((sub, i) => {
      const axisNum = i === 0 ? '' : String(i + 1);
      const { bgPts, posPts, negPts, yOf } = volcanoTracesFor(m, sub, dl, thresh);
      traces.push({ x: bgPts.map(p => p.b), y: bgPts.map(yOf), text: bgPts.map(p => p.g),
                    mode: 'markers', type: 'scattergl', name: 'Not sig', marker: { color: '#ccc', size: 4 },
                    xaxis: 'x' + axisNum, yaxis: 'y' + axisNum, showlegend: i === 0, hoverinfo: 'text+x+y' });
      traces.push({ x: negPts.map(p => p.b), y: negPts.map(yOf), text: negPts.map(p => p.g),
                    mode: 'markers', type: 'scattergl', name: dl.neg, marker: { color: dl.negColor, size: 5 },
                    xaxis: 'x' + axisNum, yaxis: 'y' + axisNum, showlegend: i === 0, hoverinfo: 'text+x+y' });
      traces.push({ x: posPts.map(p => p.b), y: posPts.map(yOf), text: posPts.map(p => p.g),
                    mode: 'markers', type: 'scattergl', name: dl.pos, marker: { color: dl.posColor, size: 5 },
                    xaxis: 'x' + axisNum, yaxis: 'y' + axisNum, showlegend: i === 0, hoverinfo: 'text+x+y' });
      annotations.push({ text: abbr(sub), showarrow: false, xref: 'x' + axisNum + ' domain', yref: 'y' + axisNum + ' domain',
                          x: 0.5, y: 1.18, font: { size: 12 } });
    });
    Plotly.newPlot(container, traces, {
      grid: { rows, columns: cols, pattern: 'independent', xgap: 0.35, ygap: 0.5 },
      title: abbr(parent) + ' \u2014 all subclusters', annotations, margin: { t: 60, b: 40 }, showlegend: true,
      height: sizeContainerForGrid(container, rows)
    }, { responsive: true });
    return;
  }

  container.style.height = '640px';
  const unit = unitSel.value;
  const { bgPts, posPts, negPts, yOf } = volcanoTracesFor(m, unit, dl, thresh);
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
  document.getElementById('plot-heatmap').style.height = '640px';

  const checked = checkedRegions();
  const regionIdx = m.heatmap.regions.map((r, i) => checked.includes(r) ? i : -1).filter(i => i >= 0);
  const regionsFiltered = regionIdx.map(i => m.heatmap.regions[i]);

  let units = m.heatmap.units;
  let z = metric.z.map(row => regionIdx.map(i => row[i]));
  const unit = unitSel.value;
  let titleScope = 'All units';
  let scopedUnits = units, scopedZ = z;

  if (unit.startsWith('PARENT:')) {
    const parent = unit.slice('PARENT:'.length);
    const keepIdx = units.map((u, i) => m.unit_to_parent[u] === parent ? i : -1).filter(i => i >= 0);
    scopedUnits = keepIdx.map(i => units[i]);
    scopedZ = keepIdx.map(i => z[i]);
    titleScope = abbr(parent) + ' subclusters';
  } else if (unit !== '__ALL__') {
    const keepIdx = m.heatmap.units.map((u, i) => u === unit ? i : -1).filter(i => i >= 0);
    scopedUnits = keepIdx.map(i => m.heatmap.units[i]);
    scopedZ = keepIdx.map(i => z[i]);
    titleScope = abbr(unit);
  }

  const aggMode = document.getElementById('heatmapAggSel').value;
  let finalUnits = scopedUnits, finalZ = scopedZ;
  if (aggMode !== 'none' && (unit === '__ALL__' || unit.startsWith('PARENT:'))) {
    const groupOf = unit === '__ALL__' ? (u => m.unit_to_parent[u]) : (() => (unit.startsWith('PARENT:') ? unit.slice('PARENT:'.length) : u));
    const agg = aggregateRows(scopedUnits, scopedZ, groupOf, aggMode);
    finalUnits = agg.units; finalZ = agg.z;
    titleScope += (aggMode === 'mean' ? ' (mean per cell type)' : ' (median per cell type)');
  }

  const label = metricLabel(metricKey) || metric.label;
  const flat = finalZ.flat().filter(v => v !== null && v !== undefined);
  const yLabels = finalUnits.map(abbr);
  const showValues = document.getElementById('heatmapShowValues').checked;
  const textMatrix = finalZ.map(row => row.map(fmtVal));

  let trace;
  if (metric.kind === 'sequential') {
    const maxV = flat.length ? Math.max(...flat) : 1;
    const seqColors = (metricKey === 'n_decrease')
      ? [[0, '#ffffff'], [1, '#2166ac']]
      : [[0, '#ffffff'], [1, '#b2182b']];
    trace = {
      z: finalZ, x: regionsFiltered, y: finalUnits, type: 'heatmap',
      colorscale: seqColors, zmin: 0, zmax: maxV || 1,
      colorbar: { title: label }
    };
  } else {
    const bound = flat.length ? Math.max(...flat.map(Math.abs)) : 1;
    trace = {
      z: finalZ, x: regionsFiltered, y: finalUnits, type: 'heatmap',
      colorscale: [[0, '#2166ac'], [0.5, '#ffffff'], [1, '#b2182b']],
      zmin: -bound, zmax: bound,
      colorbar: { title: label }
    };
  }
  trace.hovertemplate = '%{y} \u00d7 %{x}<br>' + label + ': %{z}<extra></extra>';
  if (showValues) {
    trace.text = textMatrix;
    trace.texttemplate = '%{text}';
    trace.textfont = { size: 9, color: '#111' };
  }

  Plotly.react('plot-heatmap', [trace], {
    title: label + ' \u2014 ' + titleScope,
    margin: { t: 40, l: 140 },
    plot_bgcolor: 'white', paper_bgcolor: 'white',
    xaxis: { showgrid: false, zeroline: false },
    yaxis: { tickvals: finalUnits, ticktext: yLabels, autorange: 'reversed', showgrid: false, zeroline: false }
  }, {responsive: true});

  document.getElementById('note-heatmap').textContent =
    metricKey.startsWith('raw_beta')
      ? 'Raw per-gene q-values are known to be poorly calibrated \u2014 interpret this metric cautiously relative to the mash-based ones.'
      : 'Blank cells = that unit x region combination doesn\'t exist in the data. Uses the selected LFSR threshold (' + thresh + ').';
}

function showUnavailable() {
  const msg = '<div class="msg">This model hasn\'t been (re)generated yet \u2014 check back once the pipeline finishes.</div>';
  document.getElementById('plot-bar').innerHTML = msg;
  document.getElementById('note-bar').textContent = '';
  document.getElementById('plot-volcano').innerHTML = msg;
  document.getElementById('plot-heatmap').innerHTML = msg;
  document.getElementById('note-heatmap').textContent = '';
  document.getElementById('gene-result').innerHTML = msg;
}

const CHECKPOINTS_DIR = '/scratch/easmit31/dispersion/dglm/disp_age__raw_counts';

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
  if (unit === '__ALL__' || unit.startsWith('PARENT:')) {
    container.innerHTML = '<div class="msg">Pick one specific subcluster (not "All" or a grouped cell type) to search a gene within it.</div>';
    return;
  }
  const byRegion = m.volcano[unit] || {};
  const checked = checkedRegions();
  const hits = [];
  Object.keys(byRegion).forEach(region => {
    if (!checked.includes(region)) return;
    byRegion[region].forEach(p => {
      if (p.g.toLowerCase() === query.toLowerCase()) hits.push({ region, b: p.b, l: p.l });
    });
  });

  container.innerHTML = '';
  if (hits.length === 0) {
    const msgDiv = document.createElement('div');
    msgDiv.className = 'msg';
    msgDiv.textContent = '"' + query + '" not found in the visualized sample for ' + abbr(unit) +
      ' in this model (within the currently checked regions). Either it\'s not significant here, or it fell outside the random background sample embedded for file-size reasons.';
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

  const ctDir = CHECKPOINTS_DIR + '/' + unit.replace(/_\d+$/, '') + '/dglm_checkpoints_cutoff0.5';
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
unitSel.addEventListener('change', () => { updateControlVisibility(); drawActive(); });
regionSel.addEventListener('change', drawActive);
document.getElementById('heatmapMetricSel').addEventListener('change', drawActive);
document.getElementById('heatmapShowValues').addEventListener('change', drawActive);
document.getElementById('heatmapAggSel').addEventListener('change', drawActive);
document.getElementById('barSourceSel').addEventListener('change', drawActive);
document.getElementById('barGroupModeSel').addEventListener('change', drawActive);
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
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    v2 = build_model('', 'master_dglm_combined.tsv',
                      'Age-only dispersion model, sqrt Shat, per-cell-type mashr, louvain resolution')
    v3 = build_model('_agesex', 'master_dglm_age_disp_combined.tsv',
                      'Age+sex dispersion model (age term), sqrt Shat, per-cell-type mashr, louvain resolution')

    global BASE_DIR
    BASE_DIR = '/scratch/easmit31/dispersion/dglm/disp_age__min300cells_min30animals'
    v_filtered300 = build_model('', 'master_dglm_combined.tsv',
                      'Age-only, raw Shat, mash_1by1 strong subset, >=300 cells/>=30 animals condition filter (GABAergic_neurons, microglia so far -- other cell types pending)')
    BASE_DIR = '/scratch/easmit31/dispersion/dglm/disp_age__min100cells_min30animals'
    v_filtered100 = build_model('', 'master_dglm_combined.tsv',
                      'Age-only, raw Shat, mash_1by1 strong subset, >=100 cells/>=30 animals condition filter (GABAergic_neurons, microglia so far -- other cell types pending)')
    BASE_DIR = '/scratch/easmit31/dispersion/dglm/disp_age__raw_counts'

    result = {'v2_sqrt_louvain': v2, 'v3_sqrt_louvain_agesex': v3,
              'v_filtered300_raw_lfsr_single': v_filtered300,
              'v_filtered100_raw_lfsr_single': v_filtered100}
    data_json_text = json.dumps(result)

    marker = "const modelSel = document.getElementById('modelSel');"
    idx = TEMPLATE_HEAD.find(marker)
    if idx == -1:
        raise SystemExit('Could not find injection point in template')
    html = TEMPLATE_HEAD[:idx] + f'const DATA = {data_json_text};\n\n' + TEMPLATE_HEAD[idx:]

    with open(args.out, 'w') as f:
        f.write(html)
    print(f'Saved: {args.out}')

if __name__ == '__main__':
    main()
