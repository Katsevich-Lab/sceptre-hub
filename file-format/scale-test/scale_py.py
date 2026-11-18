"""Scale test, Python side. Stage `convert`: the Gasperini at-scale data, held
in the October layout, written as a data file of the format; validated; read
back in memory and feature by feature. Stage `compare`: the file R wrote from
it, checked against the Python-written one."""
import json, os, sys, time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, "format")
import sceptre_format as sf

ROOT = Path(".")
SRC = Path(os.environ["SCEPTRE_SCALE_SOURCE"])   # a folder with screen.h5mu and grna_targets.parquet in the October layout
OUT = ROOT / "out"
results_file = OUT / "results-python.json"
results = json.loads(results_file.read_text()) if results_file.exists() else {}


def peak_gb():
    for line in open("/proc/self/status"):
        if line.startswith("VmHWM"):
            return round(int(line.split()[1]) / 1e6, 1)


def step(name, fn):
    try:
        open("/proc/self/clear_refs", "w").write("5")   # reset the peak
    except OSError:
        pass
    start = time.time()
    value = fn()
    results[name] = {"seconds": round(time.time() - start, 1), "peak_gb": peak_gb()}
    print(name, results[name], flush=True)
    results_file.write_text(json.dumps(results, indent=2))
    return value


def size_gb(path):
    return round(Path(path).stat().st_size / 1e9, 2)


stage = sys.argv[1]
if stage == "convert":
    import anndata as ad
    import mudata as md

    def load():
        return md.read_h5mu(SRC / "screen.h5mu")

    old = step("read the October-layout file (in memory)", load)
    obs = old["response"].obs.drop(columns=[c for c in old["response"].obs if c.startswith("log1p_")])
    parts = sf.split_covariates(obs, ("response", "grna"))
    parts["response"] = parts["response"].assign(n_nonzero=parts["response"]["n_nonzero"].astype(np.int32))
    parts["grna"] = parts["grna"].assign(n_nonzero=parts["grna"]["n_nonzero"].astype(np.int32))
    targets = pd.read_parquet(SRC / "grna_targets.parquet").set_index("grna_id").loc[list(old["grna"].var_names)]
    targets.index.name = None
    targets = targets.drop(columns="grna_group").rename(columns=sf.GRNA_COORDINATES).astype({"grna_target": object, "target_chr": object})
    new = md.MuData({
        "rna": ad.AnnData(X=old["response"].X, obs=parts["response"], var=old["response"].var),
        "grna": ad.AnnData(X=old["grna"].X, obs=parts["grna"], var=targets),
    })
    for name in parts["shared"].columns:
        new.obs[name] = parts["shared"][name].array
    new.uns["sceptre"] = {"moi": "high"}
    results["dataset"] = {
        "cells": int(new.n_obs), "responses": int(new["rna"].n_vars), "grnas": int(new["grna"].n_vars),
        "rna_nonzeros": int(new["rna"].X.nnz), "grna_nonzeros": int(new["grna"].X.nnz),
        "rna_count_sum": float(new["rna"].X.sum()), "grna_count_sum": float(new["grna"].X.sum()),
        "october_file_gb": size_gb(SRC / "screen.h5mu"),
    }
    out = OUT / "screen-python.h5mu"
    out.unlink(missing_ok=True)
    step("write the data file", lambda: sf.write_screen(new, out))
    results["dataset"]["format_file_gb"] = size_gb(out)

    # an assignment table of real size: every (gRNA, cell) with at least 5 gRNA UMIs
    def assignments():
        counts = new["grna"].X.tocoo()
        keep = counts.data >= 5
        table = pd.DataFrame({
            "grna_id": np.asarray(new["grna"].var_names, dtype=object)[counts.col[keep]],
            "cell_id": np.asarray(new["grna"].obs_names, dtype=object)[counts.row[keep]],
        })
        return sf.write_table(table, OUT / "assignments-python.parquet", ["grna_id", "cell_id"])

    rows = step("write an assignment table", assignments)
    results["assignments"] = {"rows": int(rows), "file_mb": round((OUT / "assignments-python.parquet").stat().st_size / 1e6, 1)}
    del old, new

    problems = step("validate the data file", lambda: sf.validate_screen(out))
    results["validation"] = {"level": sf.conformance_level(problems), "problems": problems.to_dict("records")}

    back = step("read the data file (in memory)", lambda: sf.read_screen(out))
    d = results["dataset"]
    results["read back"] = {
        "same shape and nonzeros": [back.n_obs, back["rna"].n_vars, back["grna"].n_vars, back["rna"].X.nnz, back["grna"].X.nnz]
        == [d["cells"], d["responses"], d["grnas"], d["rna_nonzeros"], d["grna_nonzeros"]],
        "same count sums": float(back["rna"].X.sum()) == d["rna_count_sum"] and float(back["grna"].X.sum()) == d["grna_count_sum"],
    }
    # per-feature sums, for R to compare with
    np.save(OUT / "rna-feature-sums.npy", np.asarray(back["rna"].X.sum(axis=0)).reshape(-1))
    rng = np.random.default_rng(1)
    picked = np.sort(rng.choice(back["rna"].n_vars, 200, replace=False))
    in_memory = [(back["rna"].X.indices[back["rna"].X.indptr[j]:back["rna"].X.indptr[j + 1]].copy(),
                  back["rna"].X.data[back["rna"].X.indptr[j]:back["rna"].X.indptr[j + 1]].copy()) for j in picked]
    del back
    json.dump(picked.tolist(), open(OUT / "picked-features.json", "w"))

    def by_feature():
        ok = True
        with h5py.File(out, "r") as f:
            X = f["mod/rna/X"]
            pointers = X["indptr"][()]
            for j, (cells, counts) in zip(picked, in_memory):
                lo, hi = int(pointers[j]), int(pointers[j + 1])
                ok &= np.array_equal(X["indices"][lo:hi], cells) and np.array_equal(X["data"][lo:hi], counts)
        return bool(ok)

    results["by feature"] = {"200 features read from disk equal the in-memory ones": step("read 200 features from disk, one at a time", by_feature)}
    step("read the assignment table", lambda: sf.read_table(OUT / "assignments-python.parquet"))

if stage == "compare":
    theirs, mine = OUT / "screen-r.h5mu", OUT / "screen-python.h5mu"
    problems = step("validate the R-written data file", lambda: sf.validate_screen(theirs))
    out = {"R-written file": {"level": sf.conformance_level(problems), "problems": problems.to_dict("records"), "gb": size_gb(theirs)}}

    def types(path):
        with h5py.File(path) as f:
            return {m: {k: str(f[f"mod/{m}/X/{k}"].dtype) for k in ("data", "indices", "indptr")} for m in f["mod"]}

    def same_arrays():
        ok = {}
        with h5py.File(theirs) as a, h5py.File(mine) as b:
            for m in ("rna", "grna"):
                for k in ("data", "indices", "indptr"):
                    x, y = a[f"mod/{m}/X/{k}"], b[f"mod/{m}/X/{k}"]
                    same = x.shape == y.shape
                    for lo in range(0, x.shape[0], 50_000_000):
                        same = same and np.array_equal(x[lo:lo + 50_000_000], y[lo:lo + 50_000_000])
                    ok[f"{m} {k}"] = bool(same)
        return ok

    out["same number types"] = types(theirs) == types(mine)
    out["count arrays identical"] = step("compare the count arrays of the two files", same_arrays)

    def same_tables():
        a, b = sf.read_screen(theirs), sf.read_screen(mine)
        return {
            "cell tables": bool(a.obs.equals(b.obs) and all(a[m].obs.equals(b[m].obs) for m in ("rna", "grna"))),
            "feature tables": bool(all(a[m].var.equals(b[m].var) for m in ("rna", "grna"))),
            "identifiers": bool(all(list(a[m].obs_names) == list(b[m].obs_names) and list(a[m].var_names) == list(b[m].var_names) for m in ("rna", "grna"))),
        }

    out["tables identical"] = step("read both files and compare their tables", same_tables)
    r_table = OUT / "assignments-r.parquet"
    if r_table.exists():
        out["assignment table written by R equals Python's"] = bool(sf.same_table(sf.read_table(r_table), sf.read_table(OUT / "assignments-python.parquet"))) if False else \
            bool(sf.read_table(r_table).equals(sf.read_table(OUT / "assignments-python.parquet")))
    results["compare"] = out
    print(json.dumps(out, indent=2))
results_file.write_text(json.dumps(results, indent=2))
