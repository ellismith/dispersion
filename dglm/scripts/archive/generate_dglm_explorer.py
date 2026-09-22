#!/usr/bin/env python3
"""
Patch the existing DGLM explorer generator so its louvain candidates discover
all per-parent-cell-type master TSVs automatically.

Run:
    python generate_dglm_explorer.py
"""

import os
import re
from pathlib import Path

SCRIPT = Path("/scratch/easmit31/dispersion/dglm/scripts/generate_dglm_explorer.py")

text = SCRIPT.read_text()

raw_root = (
    "/scratch/easmit31/dispersion/dglm/"
    "disp_age__mean_full_louvain_percell_raw"
)
sqrt_root = (
    "/scratch/easmit31/dispersion/dglm/"
    "disp_age__mean_full_louvain_percell_sqrt"
)

new_block = f"""LOUVAIN_RAW_ROOT = (
    '{raw_root}'
)
LOUVAIN_SQRT_ROOT = (
    '{sqrt_root}'
)


def louvain_master_paths(root):
    \"\"\"Discover every per-parent-cell-type louvain master TSV on disk.\"\"\"
    if not os.path.isdir(root):
        return []

    paths = []
    for entry in sorted(os.scandir(root), key=lambda x: x.name.lower()):
        if not entry.is_dir():
            continue
        master = os.path.join(entry.path, 'master_dglm_combined.tsv')
        if os.path.isfile(master) and os.path.getsize(master) > 0:
            paths.append(master)

    return paths
"""

text = re.sub(
    r"LOUVAIN_CELL_TYPES\s*=\s*\[.*?\]\s*\n\s*LOUVAIN_RAW_ROOT\s*=\s*\(.*?\)\s*\n\s*LOUVAIN_SQRT_ROOT\s*=\s*\(.*?\)\s*",
    new_block + "\n",
    text,
    flags=re.DOTALL,
)

text = re.sub(
    r"def louvain_master_paths\(root\):\n.*?(?=\n\ndef load_candidate)",
    "",
    text,
    flags=re.DOTALL,
)

insert_at = text.find("\ndef load_candidate")
if insert_at == -1:
    raise RuntimeError("Could not find load_candidate() in generate_dglm_explorer.py")

if "def louvain_master_paths(root):" not in text:
    text = text[:insert_at] + "\n\n" + new_block + text[insert_at:]

text = text.replace(
    "disp_age__mean_full_louvain_percell_SQRT",
    "disp_age__mean_full_louvain_percell_sqrt",
)

text = text.replace(
    "v2_sqrt_louvain: '/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain_percell_SQRT'",
    "v2_sqrt_louvain: '/scratch/easmit31/dispersion/dglm/disp_age__mean_full_louvain_percell_sqrt'",
)

SCRIPT.write_text(text)

print(f"Patched: {SCRIPT}")
print(f"Raw root:  {raw_root}")
print(f"Sqrt root: {sqrt_root}")
