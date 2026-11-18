"""Readers we did not write. The example files are opened with tools outside the
prototype (anndata and scanpy for a modality; polars and duckdb for the tables)
and what they return is compared with what the files say and with pyarrow.

    python verify/independent_readers.py
"""
import sys
import warnings
from pathlib import Path

import anndata as ad
import duckdb
import h5py
import numpy as np
import polars as pl
import pyarrow.parquet as pq
import scanpy as sc

sys.path.insert(0, "verify")
from common import Results

warnings.simplefilter("ignore")
results = Results("independent-readers")
FILES = {
    "example/highmoi/screen.h5mu": "written by R",
    "example/lowmoi/screen.h5mu": "written by R",
    "conformance/valid/highmoi_written_by_python.h5mu": "written by Python",
}
for file, writer in FILES.items():
    with h5py.File(file) as f:
        for modality in ("rna", "grna"):
            # anndata reads one modality as an AnnData, with no MuData and none of our code
            adata = ad.io.read_elem(f[f"mod/{modality}"])
            shape = tuple(int(v) for v in f[f"mod/{modality}/X"].attrs["shape"])
            results.add(f"anndata reads {modality} of {file} ({writer})", adata.shape == shape, f"{adata.shape}")
            # scanpy on the matrix exactly as anndata returns it, with nothing done to it:
            # both index arrays are stored in one width (S4), which scipy requires in some
            # operations. percent_top=None because the default asks for the top 500
            # features, which a 100-feature example does not have.
            results.add(f"both index arrays of {modality} of {file} arrive in one width", adata.X.indices.dtype == adata.X.indptr.dtype,
                        f"{adata.X.indices.dtype}, {adata.X.indptr.dtype}")
            adata.X.sum(); adata[:10, :5].X.toarray(); adata.X.tocsr()
            metrics = sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False)[0]
            results.add(f"scanpy's total counts equal the stored n_umis, {modality} of {file}",
                        np.array_equal(metrics["total_counts"].to_numpy(), adata.obs["n_umis"].to_numpy(dtype=float)))
            results.add(f"scanpy's number of detected features equals the stored n_nonzero, {modality} of {file}",
                        np.array_equal(metrics["n_genes_by_counts"].to_numpy(), adata.obs["n_nonzero"].to_numpy()))

for folder in ("example/highmoi/analyses/left_sided", "example/highmoi/analyses/singleton", "example/lowmoi/analyses/default",
               "conformance/valid/analyses/left_sided_written_by_python"):
    for path in sorted(Path(folder).rglob("*.parquet")):
        if "fits" in path.parts:
            continue
        name = f"{path.relative_to(Path(folder).parent)}"
        arrow = pq.read_table(path)
        expected = arrow.to_pylist()
        polars = pl.read_parquet(path)
        results.add(f"polars reads {name} to the same rows as pyarrow", polars.columns == arrow.column_names and polars.to_dicts() == expected)
        duck = duckdb.sql(f"select * from read_parquet('{path}')").to_arrow_table()
        results.add(f"duckdb reads {name} to the same rows as pyarrow", duck.column_names == arrow.column_names and duck.to_pylist() == expected)
    # the fit cache is a folder of part files: both read it as one table
    files = sorted(Path(folder).glob("fits/response/*.parquet"))
    if files:
        parts = f"{folder}/fits/response/*.parquet"
        rows = sum(pq.ParquetFile(p).metadata.num_rows for p in files)
        results.add(f"polars reads the part files of the fit cache of {folder} as one table", pl.read_parquet(parts).height == rows)
        results.add(f"duckdb reads the part files of the fit cache of {folder} as one table",
                    duckdb.sql(f"select count(*) from read_parquet('{parts}')").fetchall()[0][0] == rows)
results.save()
