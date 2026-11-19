"""The 64-bit branch of D4. No test file is large enough to need 64-bit index
arrays, so the writer's limit is lowered and the highmoi example is written with
both index arrays 64-bit, as a file past the real limit would have them. Python
then reads it back with anndata, scipy and scanpy; verify/wide_index.R reads it
with R's readers. This shows the types are handled. It cannot show that values
above 2,147,483,647 are, which needs a file of that size.

    python verify/wide_index.py
"""

import sys
from pathlib import Path

import anndata as ad
import h5py
import mudata as md
import numpy as np
import scanpy as sc

sys.path.insert(0, "format")
sys.path.insert(0, "verify")
import perturb_seq_data as sf
from common import Results

results = Results("wide-index")
source = "example/highmoi/screen.h5mu"
out = Path("verify/out/wide")
out.mkdir(parents=True, exist_ok=True)
wide = out / "screen-wide.h5mu"
wide.unlink(missing_ok=True)
sf.write_screen(md.read_h5mu(source), wide, index_limit=1)
with h5py.File(wide) as f, h5py.File(source) as g:
    for modality in ("rna", "grna"):
        types = {a: str(f[f"mod/{modality}/X/{a}"].dtype) for a in ("data", "indices", "indptr")}
        results.add(f"written with the limit lowered, {modality} has 64-bit index arrays and 32-bit counts",
                    types == {"data": "int32", "indices": "int64", "indptr": "int64"}, f"{types}")
        adata, narrow = ad.io.read_elem(f[f"mod/{modality}"]), ad.io.read_elem(g[f"mod/{modality}"])
        results.add(f"anndata reads {modality} with both index arrays 64-bit", (adata.X.indices.dtype, adata.X.indptr.dtype) == (np.int64, np.int64))
        results.add(f"the counts of {modality} equal the 32-bit file's", (adata.X != narrow.X).nnz == 0)
        metrics = sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False)[0]
        results.add(f"scanpy's total counts equal the stored n_umis, {modality}",
                    np.array_equal(metrics["total_counts"].to_numpy(), adata.obs["n_umis"].to_numpy(dtype=float)))
results.add("the validator passes it under the lowered limit", len(sf.validate_screen(wide, index_limit=1)) == 0)
found = sf.validate_screen(wide)
results.add("under the real limit it does not conform, breaking D4.2 and D4.3 only", sorted(set(found["rule"])) == ["D4.2", "D4.3"])
back = md.read_h5mu(wide)
results.add("mudata returns the same counts", all((back[m].X != md.read_h5mu(source)[m].X).nnz == 0 for m in ("rna", "grna")))
results.save()
