"""Random inputs, seeded and small: a dataset and an analysis in a neutral form
(JSON), which R and Python each write as files of the format. Run first; it also
does the Python writing.

    python verify/random_inputs.py
"""
import json
import shutil
import sys

import numpy as np

sys.path.insert(0, "verify")
sys.path.insert(0, "format")
import sceptre_format as sf
from common import OUT
from two_writers import write_all

SEEDS = (1, 2, 3)
ROOT = OUT / "random"
shutil.rmtree(ROOT, ignore_errors=True)
for seed in SEEDS:
    rng = np.random.default_rng(seed)
    n_cells, n_genes, n_targets = int(rng.integers(30, 80)), int(rng.integers(5, 15)), int(rng.integers(2, 5))
    cells = [f"cell_{i}" for i in range(n_cells)]
    genes = [f"gène_{i}" for i in range(n_genes)]   # identifiers outside ASCII
    targets = [f"site_{i}" for i in range(n_targets)]
    grnas = [f"{t}_g{j}" for t in targets for j in (1, 2)] + ["nt_1", "nt_2"]
    grna_targets = [t for t in targets for _ in (1, 2)] + ["non-targeting"] * 2

    def triplets(n_features, density, high):
        mask = rng.random((n_cells, n_features)) < density
        rows, cols = np.nonzero(mask)
        return [[int(r), int(c), int(v)] for r, c, v in zip(rows, cols, rng.integers(1, high, len(rows)))]

    grna_counts = triplets(len(grnas), 0.25, 40)
    assignments = sorted({(grnas[c], cells[r]) for r, c, v in grna_counts if v >= 5})
    in_use = [bool(v) for v in rng.random(n_cells) > 0.1]
    pairs = sorted({(t, g) for t in targets for g in genes if rng.random() < 0.6})
    p_values = [None if rng.random() < 0.1 else float(rng.random()) for _ in pairs]
    inputs = {
        "seed": seed, "id": f"random-{seed}", "moi": "high", "cells": cells, "genes": genes, "grnas": grnas,
        "rna": triplets(n_genes, 0.3, 200), "grna": grna_counts,
        "batch": {"levels": ["b2", "b1"], "values": [["b2", "b1"][int(v)] for v in rng.integers(0, 2, n_cells)]},
        "quality": [None if rng.random() < 0.1 else float(rng.normal()) for _ in cells],
        "grna_target": grna_targets,
        "target_chr": [f"chr{1 + i % 3}" for i in range(2 * n_targets)] + [None, None],
        "target_start": [1000 * (i + 1) for i in range(2 * n_targets)] + [None, None],
        "target_end": [1000 * (i + 1) + 19 for i in range(2 * n_targets)] + [None, None],
        "assignments": [list(a) for a in assignments], "in_use": in_use,
        "pairs": [list(p) for p in pairs],
        "results": [{"grna_target": t, "response_id": g, "pass_qc": p is not None, "p_value": p,
                     "fold_change": None if p is None else float(rng.lognormal()), "se_fold_change": None if p is None else float(rng.random())}
                    for (t, g), p in zip(pairs, p_values)],
        "fits": [{"response_id": g, "coefficients": [float(v) for v in rng.normal(size=3)], "theta": float(rng.lognormal())} for g in genes],
    }
    folder = ROOT / f"seed-{seed}"
    folder.mkdir(parents=True)
    (folder / "input.json").write_text(json.dumps(inputs))
    write_all(inputs, folder, "python")
print(f"random inputs for seeds {SEEDS}; Python wrote its files")
