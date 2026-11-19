"""Python steps of the past-the-limit test. Run as
    python py_steps.py {write|validate|read|disk|tools} ROOT
ROOT holds code/ (format/, spec/), data/ and results/."""

import json
import resource
import sys
import threading
import time
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np

step, ROOT = sys.argv[1], Path(sys.argv[2])
sys.path.insert(0, str(ROOT / "code" / "format"))
sys.path.insert(0, str(ROOT / "code" / "scripts"))
import gen  # noqa: E402

FILE = ROOT / "data" / "screen.h5mu"
TRUTH = ROOT / "data" / "truth.npz"
SAMPLES = ROOT / "data" / "truth_samples.h5"
OUT = ROOT / "results" / f"py-{step}.json"
LIMIT = 2**31 - 1
record = {"step": step, "file": str(FILE), "checks": [], "timings": {}, "started": datetime.now().isoformat()}


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def mem():
    status = dict(l.split(":", 1) for l in Path("/proc/self/status").read_text().splitlines() if ":" in l)
    return status["VmRSS"].strip(), status["VmHWM"].strip()


def log(*a):
    print(f"[{now()}]", *a, flush=True)


def heartbeat():
    while True:
        time.sleep(60)
        rss, hwm = mem()
        log(f"heartbeat rss={rss} peak={hwm}")


threading.Thread(target=heartbeat, daemon=True).start()


def save():
    record["peak_rss_kb_getrusage"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    record["finished"] = datetime.now().isoformat()
    OUT.write_text(json.dumps(record, indent=1, default=str))


def check(name, ok, detail=""):
    record["checks"].append({"check": name, "ok": bool(ok), "detail": str(detail)})
    log("PASS" if ok else "FAIL", name, detail)
    save()


class timed:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t = time.time()
        log("start", self.name)

    def __exit__(self, kind, value, tb):
        dt = time.time() - self.t
        record["timings"][self.name] = {"seconds": round(dt, 1), "ok": kind is None, "peak_rss_after": mem()[1]}
        if kind is not None:
            record["timings"][self.name]["error"] = "".join(traceback.format_exception(kind, value, tb))
            log("ERROR in", self.name, f"after {dt:.1f}s:\n" + record["timings"][self.name]["error"])
        else:
            log(f"done {self.name} in {dt:.1f}s, peak rss {mem()[1]}")
        save()
        return True   # record and continue


def sample_genes(indptr):
    """first, last, the genes whose stored values cross positions 2^31-1, and random."""
    straddle = int(np.searchsorted(indptr, LIMIT, side="right") - 1)
    rng = np.random.default_rng(7)
    picks = [0, 1, 2, straddle - 1, straddle, straddle + 1, gen.N_GENES - 3, gen.N_GENES - 2, gen.N_GENES - 1]
    picks += sorted(rng.choice(gen.N_GENES, 11, replace=False).tolist())
    return sorted({p for p in picks if 0 <= p < gen.N_GENES}), straddle


def per_gene(indptr, data):
    """per-gene sums (int64) of a csc matrix's arrays, empty genes included."""
    sums = np.zeros(len(indptr) - 1, dtype=np.int64)
    for g in range(len(indptr) - 1):
        a, b = int(indptr[g]), int(indptr[g + 1])
        sums[g] = data[a:b].sum(dtype=np.int64)
    return sums


def compare_truth(label, indptr, data, truth):
    nnz = np.diff(np.asarray(indptr, dtype=np.int64))
    check(f"{label}: nnz per gene equals ground truth", np.array_equal(nnz, truth["gene_nnz"]),
          f"total nnz {int(nnz.sum())} vs truth {int(truth['gene_nnz'].sum())}")
    sums = per_gene(indptr, data)
    check(f"{label}: sum per gene equals ground truth", np.array_equal(sums, truth["gene_sum"]),
          f"genes differing: {int((sums != truth['gene_sum']).sum())}; total {int(sums.sum())} vs {int(truth['gene_sum'].sum())}")


# ------------------------------------------------------------------ write
if step == "write":
    import h5py
    import mudata as md
    import anndata as ad
    import pandas as pd
    from scipy import sparse
    import perturb_seq_data as sf

    with timed("generate"):
        nnz = np.array([gen.gene_nnz(g) for g in range(gen.N_GENES)], dtype=np.int64)
        indptr = np.concatenate([[0], np.cumsum(nnz)]).astype(np.int64)
        total = int(indptr[-1])
        log(f"total nnz {total} (limit {LIMIT}, over by {total - LIMIT})")
        indices = np.empty(total, dtype=np.int64)
        data = np.empty(total, dtype=np.int32)
        gene_sum = np.zeros(gen.N_GENES, dtype=np.int64)
        cell_sum = np.zeros(gen.N_CELLS, dtype=np.int64)
        cell_nnz = np.zeros(gen.N_CELLS, dtype=np.int64)
        for g in range(gen.N_GENES):
            pos, vals = gen.gene(g)
            assert len(pos) == nnz[g]
            a, b = indptr[g], indptr[g + 1]
            indices[a:b] = pos
            data[a:b] = vals
            gene_sum[g] = vals.sum(dtype=np.int64)
            cell_sum[pos] += vals
            cell_nnz[pos] += 1
            if g % 500 == 0:
                log(f"generated gene {g}/{gen.N_GENES}, rss {mem()[0]}")
        np.savez(TRUTH, gene_nnz=nnz, gene_sum=gene_sum, cell_sum=cell_sum, cell_nnz=cell_nnz)
        picks, straddle = sample_genes(indptr)
        with h5py.File(SAMPLES, "w") as f:
            f["genes"] = np.array(picks, dtype=np.int64)
            for g in picks:
                pos, vals = gen.gene(g)
                f[f"pos/{g}"] = pos
                f[f"val/{g}"] = vals
        import csv
        with open(ROOT / "data" / "truth_genes.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["gene", "gene_id", "nnz", "sum"])
            for g, gid in enumerate(gen.gene_ids()):
                w.writerow([g, gid, int(nnz[g]), int(gene_sum[g])])
        record["dataset"] = {"n_cells": gen.N_CELLS, "n_genes": gen.N_GENES, "rna_nnz": total, "rna_nnz_minus_limit": total - LIMIT,
                             "rna_count_total": int(gene_sum.sum()), "max_count": int(data.max()), "sample_genes": picks,
                             "straddle_gene": straddle, "straddle_range": [int(indptr[straddle]), int(indptr[straddle + 1])] if straddle < gen.N_GENES else None}
        save()
    with timed("build MuData"):
        rna = sparse.csc_matrix((data, indices, indptr), shape=(gen.N_CELLS, gen.N_GENES))
        record["dataset"]["scipy_index_dtypes_in_memory"] = [str(rna.indices.dtype), str(rna.indptr.dtype)]
        del data, indices
        cells = gen.cell_ids()
        gx, gvar = gen.grna()
        rna_ad = ad.AnnData(X=rna, obs=pd.DataFrame(index=cells), var=pd.DataFrame({"gene_name": [f"GENE{g}" for g in range(gen.N_GENES)]}, index=gen.gene_ids()))
        grna_ad = ad.AnnData(X=gx, obs=pd.DataFrame(index=cells), var=gvar)
        mdata = md.MuData({"rna": rna_ad, "grna": grna_ad})
        mdata.obs["batch"] = pd.Categorical(np.array(["lane1", "lane2", "lane3", "lane4"])[np.arange(gen.N_CELLS) % 4])
        mdata.uns["title"] = "Synthetic screen past 2^31-1 nonzero counts"
        mdata.uns["perturbation_type"] = "CRISPRi"
        mdata.uns["organism"] = "Homo sapiens"
        record["dataset"]["grna_shape"] = list(gx.shape)
        record["dataset"]["grna_nnz"] = int(gx.nnz)
        log("rna X id same as input:", mdata["rna"].X is rna, type(mdata["rna"].X))
    FILE.unlink(missing_ok=True)
    with timed("write_screen (real limit)"):
        sf.write_screen(mdata, FILE)
    del mdata, rna_ad, rna
    record["file_size_bytes"] = FILE.stat().st_size if FILE.exists() else None
    with h5py.File(FILE) as f:
        types = {m: {a: str(f[f"mod/{m}/X/{a}"].dtype) for a in ("data", "indices", "indptr")} for m in ("rna", "grna")}
        record["stored_types"] = types
        check("rna: data int32, indices int64, indptr int64", types["rna"] == {"data": "int32", "indices": "int64", "indptr": "int64"}, types["rna"])
        check("grna (small, own size under the limit): data int32, indices int32, indptr int32 per D4 (index type is per modality)",
              types["grna"] == {"data": "int32", "indices": "int32", "indptr": "int32"}, types["grna"])
        ip = f["mod/rna/X/indptr"]
        check("rna indptr ends at nnz > 2^31-1", int(ip[-1]) == record["dataset"]["rna_nnz"] and int(ip[-1]) > LIMIT, int(ip[-1]))
        check("rna shape attr", list(f["mod/rna/X"].attrs["shape"]) == [gen.N_CELLS, gen.N_GENES], list(f["mod/rna/X"].attrs["shape"]))
        record["filters"] = {a: f[f"mod/rna/X/{a}"].compression for a in ("data", "indices", "indptr")}
    save()

# ------------------------------------------------------------------ validate
elif step == "validate":
    import perturb_seq_data as sf

    with timed("validate_screen (real limit)"):
        found = sf.validate_screen(FILE)
        record["problems"] = found.to_dict("records")
        check("validate_screen reports nothing", len(found) == 0, found.to_string())

# ------------------------------------------------------------------ disk (h5py, one gene at a time)
elif step == "disk":
    import h5py

    truth = np.load(TRUTH)
    with h5py.File(SAMPLES) as s:
        picks = s["genes"][()].tolist()
    with timed("h5py: read indptr + sampled genes from disk"):
        with h5py.File(FILE) as f, h5py.File(SAMPLES) as s:
            X = f["mod/rna/X"]
            indptr = X["indptr"][()]
            record["indptr_dtype"] = str(indptr.dtype)
            per = {}
            for g in picks:
                t = time.time()
                a, b = int(indptr[g]), int(indptr[g + 1])
                pos, vals = X["indices"][a:b], X["data"][a:b]
                per[g] = round(time.time() - t, 3)
                ok = np.array_equal(pos, s[f"pos/{g}"][()]) and np.array_equal(vals, s[f"val/{g}"][()])
                check(f"h5py gene {g} (stored values {a}..{b}) equals regenerated truth exactly", ok,
                      f"nnz {b - a}, sum {int(vals.sum(dtype=np.int64))} vs {int(truth['gene_sum'][g])}")
            record["per_gene_seconds"] = per
            check("h5py indptr diff equals truth nnz", np.array_equal(np.diff(indptr), truth["gene_nnz"]))

# ------------------------------------------------------------------ read (mudata in memory)
elif step == "read":
    import mudata as md

    truth = np.load(TRUTH)
    with timed("mudata.read_h5mu (in memory)"):
        mdata = md.read_h5mu(FILE)
        X = mdata["rna"].X
        record["in_memory"] = {"type": type(X).__name__, "shape": list(X.shape), "nnz": int(X.nnz),
                               "dtypes": [str(X.data.dtype), str(X.indices.dtype), str(X.indptr.dtype)],
                               "grna_dtypes": [str(mdata["grna"].X.data.dtype), str(mdata["grna"].X.indices.dtype), str(mdata["grna"].X.indptr.dtype)]}
        log(record["in_memory"])
    with timed("compare in-memory read to ground truth"):
        compare_truth("mudata in memory", X.indptr, X.data, truth)
        cell_nnz = np.bincount(X.indices, minlength=X.shape[0])
        check("mudata in memory: nnz per cell equals ground truth", np.array_equal(cell_nnz, truth["cell_nnz"]))
        obs = mdata["rna"].obs
        check("n_umis written by write_screen equals ground-truth cell sums", np.array_equal(obs["n_umis"].to_numpy(), truth["cell_sum"].astype(float)))
        check("n_nonzero written by write_screen equals ground-truth cell nnz", np.array_equal(obs["n_nonzero"].to_numpy(), truth["cell_nnz"]))
        import h5py
        with h5py.File(SAMPLES) as s:
            for g in s["genes"][()].tolist():
                a, b = int(X.indptr[g]), int(X.indptr[g + 1])
                check(f"in-memory gene {g} equals regenerated truth exactly",
                      np.array_equal(X.indices[a:b], s[f"pos/{g}"][()]) and np.array_equal(X.data[a:b], s[f"val/{g}"][()]))

# ------------------------------------------------------------------ tools that know nothing of the format
elif step == "tools":
    import anndata as ad
    import h5py
    import scanpy as sc

    truth = np.load(TRUTH)
    record["versions"] = {"anndata": ad.__version__, "scanpy": sc.__version__}
    with timed("anndata.io.sparse_dataset: sampled genes from disk"):
        with h5py.File(FILE) as f, h5py.File(SAMPLES) as s:
            ds = ad.io.sparse_dataset(f["mod/rna/X"])
            record["sparse_dataset"] = {"format": ds.format, "shape": list(ds.shape)}
            for g in [int(x) for x in s["genes"][()]]:
                col = ds[:, g]
                ok = np.array_equal(col.indices.astype(np.int64), s[f"pos/{g}"][()]) and np.array_equal(col.data, s[f"val/{g}"][()])
                check(f"anndata sparse_dataset gene {g} equals truth", ok, f"nnz {col.nnz}")
    with timed("anndata.io.read_elem(mod/rna) in memory"):
        with h5py.File(FILE) as f:
            adata = ad.io.read_elem(f["mod/rna"])
        X = adata.X
        record["anndata"] = {"type": type(X).__name__, "nnz": int(X.nnz), "dtypes": [str(X.data.dtype), str(X.indices.dtype), str(X.indptr.dtype)]}
        log(record["anndata"])
    with timed("compare anndata read to truth"):
        compare_truth("anndata", X.indptr, X.data, truth)
    with timed("scanpy calculate_qc_metrics"):
        obs_m, var_m = sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=False)
        check("scanpy var total_counts equals truth gene sums", np.array_equal(var_m["total_counts"].to_numpy(), truth["gene_sum"].astype(float)),
              f"max abs diff {np.abs(var_m['total_counts'].to_numpy() - truth['gene_sum']).max()}")
        check("scanpy var n_cells_by_counts equals truth gene nnz", np.array_equal(var_m["n_cells_by_counts"].to_numpy(), truth["gene_nnz"]))
        check("scanpy obs total_counts equals truth cell sums", np.array_equal(obs_m["total_counts"].to_numpy(), truth["cell_sum"].astype(float)),
              f"max abs diff {np.abs(obs_m['total_counts'].to_numpy() - truth['cell_sum']).max()}")
        check("scanpy obs n_genes_by_counts equals truth cell nnz", np.array_equal(obs_m["n_genes_by_counts"].to_numpy(), truth["cell_nnz"]))
save()
log("STEP-END", step)
