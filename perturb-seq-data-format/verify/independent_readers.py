"""Readers we did not write. The example data files are opened with anndata and
scanpy, one modality at a time, with no MuData and none of this folder's code,
and what they return is compared with what the files say.

    python verify/independent_readers.py
"""
import sys
import warnings

import anndata as ad
import h5py
import numpy as np
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
            # both index arrays are stored in one width (D4), which scipy requires in some
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

results.save()
