"""Deterministic synthetic Perturb-seq counts, one gene at a time.

Every gene's cell positions and counts come from its own seeded generator, so any
gene can be regenerated exactly later without the whole matrix (the ground truth
for spot checks)."""

import numpy as np

import os

if os.environ.get("PASTLIMIT_SMOKE"):
    N_CELLS, N_GENES, TARGET_NNZ = 20_000, 500, 2_500_000
else:
    N_CELLS, N_GENES, TARGET_NNZ = 1_100_000, 8_000, 2_200_000_000
SEED = 20261119


def gene_params():
    rng = np.random.default_rng([SEED, 0])
    p = rng.beta(1.2, 3.6, N_GENES)
    p = np.clip(p * TARGET_NNZ / (N_CELLS * p.sum()), 0.0005, 0.97)
    lam = 4.0 * p / (1.0 - p + 0.05)    # denser genes have higher counts
    return p, lam


P, LAM = gene_params()


def gene_nnz(g):
    return int(np.random.default_rng([SEED, 1, g]).binomial(N_CELLS, P[g]))


def gene(g):
    """(cell positions int64, sorted; counts int32 >= 1) of gene g."""
    rng = np.random.default_rng([SEED, 1, g])
    k = int(rng.binomial(N_CELLS, P[g]))
    pos = np.sort(rng.choice(N_CELLS, k, replace=False, shuffle=False)).astype(np.int64)
    vals = (1 + rng.poisson(LAM[g], k)).astype(np.int32)
    return pos, vals


def cell_ids():
    return np.array([f"cell{i:07d}" for i in range(N_CELLS)], dtype=object)


def gene_ids():
    return np.array([f"ENSG{g:011d}" for g in range(N_GENES)], dtype=object)


def grna(n_cells=N_CELLS):
    """A small gRNA modality: 300 targeting gRNAs (3 per gene for genes 0..99)
    and 100 non-targeting. Each cell has one main gRNA, some a second one."""
    from scipy import sparse
    import pandas as pd

    rng = np.random.default_rng([SEED, 2])
    ids = [f"g{g:011d}_{j}" for g in range(100) for j in range(3)] + [f"NT_{j:03d}" for j in range(100)]
    targets = [f"ENSG{g:011d}" for g in range(100) for j in range(3)] + ["non-targeting"] * 100
    kinds = ["targeting"] * 300 + ["non-targeting"] * 100
    n = len(ids)
    main = rng.integers(0, n, n_cells)
    main_count = 5 + rng.poisson(20, n_cells)
    second_cells = np.flatnonzero(rng.random(n_cells) < 0.3)
    second = (main[second_cells] + 1 + rng.integers(0, n - 1, len(second_cells))) % n
    rows = np.concatenate([np.arange(n_cells), second_cells])
    cols = np.concatenate([main, second])
    vals = np.concatenate([main_count, 1 + rng.poisson(1, len(second_cells))]).astype(np.int32)
    x = sparse.csc_matrix((vals, (rows, cols)), shape=(n_cells, n))
    var = pd.DataFrame({"grna_target": targets, "grna_type": kinds}, index=ids)
    return x, var
