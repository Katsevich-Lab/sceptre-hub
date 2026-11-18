"""Two writers, one file. Two implementations each wrote the same random content
(verify/random_inputs.py, verify/two_writers.R). This compares the two sets of
files array by array and table by table, checks both with the validator, and
reads the first one's files back to compare them with the source, not with a
rewrite.

    python verify/two_writers.py            # the R library against the Python one
    python verify/two_writers.py A B        # any two implementations

An implementation is tested by having it read verify/out/random/seed-*/input.json
and write `screen-<name>.h5mu` and `analysis-<name>/` beside it; nothing here
imports the implementation. See verify/README.md.
"""
import json
import sys
from pathlib import Path

import h5py
import mudata as md
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, "verify")
sys.path.insert(0, "format")
import sceptre_format as sf
from common import OUT, Results

FORMULA = "~log(rna_n_umis) + batch"
TERMS = ["(Intercept)", "log(rna_n_umis)", "batchb1"]
SETTINGS = {"formula": FORMULA, "side": "both", "resampling_mechanism": "permutations", "resampling_approximation": "spa",
            "control_group": "complement", "grna_integration_strategy": "union", "B1": 499, "B2": 4999, "B3": 0,
            "multiple_testing_method": "BH", "multiple_testing_alpha": 0.1, "seed": 4, "response_fit_method": "batched", "grna_fit_method": "pysceptre"}
STAMP = "2026-01-01T00:00:00Z"


def write_all(inputs, folder, language):
    """Write the dataset and the analysis of `inputs` with the Python library."""
    import anndata as ad
    from scipy import sparse

    def counts(triplets, n_features):
        r, c, v = (np.array(x) for x in zip(*triplets))
        return sparse.csc_matrix((v, (r, c)), shape=(len(inputs["cells"]), n_features))

    cells = inputs["cells"]
    grna_var = pd.DataFrame({
        "grna_target": inputs["grna_target"], "target_chr": inputs["target_chr"],
        "target_start": pd.array(inputs["target_start"], dtype="Int32"), "target_end": pd.array(inputs["target_end"], dtype="Int32"),
    }, index=inputs["grnas"])
    mdata = md.MuData({
        "rna": ad.AnnData(X=counts(inputs["rna"], len(inputs["genes"])), obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=inputs["genes"])),
        "grna": ad.AnnData(X=counts(inputs["grna"], len(inputs["grnas"])), obs=pd.DataFrame(index=cells), var=grna_var),
    })
    mdata.obs["batch"] = pd.Categorical(inputs["batch"]["values"], categories=inputs["batch"]["levels"])
    mdata.obs["quality"] = [np.nan if v is None else v for v in inputs["quality"]]
    mdata.uns["sceptre"] = {"moi": inputs["moi"], "id": inputs["id"]}
    screen = folder / f"screen-{language}.h5mu"
    sf.write_screen(mdata, screen)

    analysis = folder / f"analysis-{language}"
    provenance = {"package": "verify", "version": "0", "language": "Python", "written_at": STAMP}
    back = md.read_h5mu(screen)
    sf.init_analysis(analysis, sf.screen_link(back), "rna")
    pairs = pd.DataFrame(inputs["pairs"], columns=["grna_target", "response_id"]).astype(object)
    none = pairs.iloc[:0]
    results = pd.DataFrame(inputs["results"])
    sf.write_stage(analysis, "parameters", {"settings": SETTINGS, "design": {"terms": TERMS, "sha256": "ab" * 32}},
                   {"parameters/discovery_pairs": pairs, "parameters/positive_control_pairs": none}, provenance)
    sf.write_stage(analysis, "assignment", {"settings": {"method": "thresholding", "hyperparameters": {"threshold": 5}}},
                   {"assignment/assignments": pd.DataFrame(inputs["assignments"], columns=["grna_id", "cell_id"]).astype(object)}, provenance)
    sf.write_stage(analysis, "qc", {"settings": {"response_n_umis_range": [0.01, 0.99], "response_n_nonzero_range": [0.01, 0.99],
                                                 "p_mito_threshold": 0.2, "pairwise_qc": True, "B3": 0}},
                   {"qc/cells": pd.DataFrame({"cell_id": np.array(cells, dtype=object), "in_use": inputs["in_use"]}),
                    "qc/discovery_pair_qc": results[["grna_target", "response_id", "pass_qc"]],
                    "qc/positive_control_pair_qc": results[["grna_target", "response_id", "pass_qc"]].iloc[:0]}, provenance)
    called = [r["p_value"] for r in inputs["results"] if r["p_value"] is not None and r["p_value"] < 0.2]
    sf.write_stage(analysis, "discovery", {"settings": {}, "significance_threshold": max(called) if called else None},
                   {"discovery/results": results.astype({"p_value": float, "fold_change": float, "se_fold_change": float})}, provenance)
    in_use = [c for c, keep in zip(cells, inputs["in_use"]) if keep]
    sf.write_fits(analysis, pd.DataFrame(inputs["fits"]), {"screen_id": inputs["id"], "formula": FORMULA, "cells_in_use_sha256": sf.sha256_ids(in_use),
                                                           "response_fit_method": "batched", "control_group": "complement"})


def column_values(node):
    """A column of a cell or feature table, decoded: its encoding and its values."""
    encoding = node.attrs.get("encoding-type")
    if encoding == "categorical":
        levels = [v.decode() if isinstance(v, bytes) else v for v in node["categories"][()]]
        return encoding, levels, [None if c < 0 else levels[c] for c in node["codes"][()]]
    values = node[()]
    if values.dtype.kind in "OS":
        return encoding, None, [v.decode() if isinstance(v, bytes) else v for v in values]
    return f"{encoding}:{values.dtype}", None, [None if isinstance(v, float) and np.isnan(v) else v for v in values.tolist()]


def same_tables(a, b):
    names = sorted(k for k in a if k != "__categories")
    if names != sorted(k for k in b if k != "__categories"):
        return f"columns differ: {names} / {sorted(b)}"
    order = lambda g: [v.decode() if isinstance(v, bytes) else v for v in g.attrs["column-order"]]
    if order(a) != order(b):
        return f"column order differs: {order(a)} / {order(b)}"
    for name in names:
        if column_values(a[name]) != column_values(b[name]):
            return f"`{name}` differs: {column_values(a[name])[0]} / {column_values(b[name])[0]}"
    return True


if __name__ == "__main__":
    results = Results("two-writers")
    A, B = (sys.argv[1:3] if len(sys.argv) >= 3 else ["r", "python"])
    for folder in sorted((OUT / "random").glob("seed-*")):
        seed = folder.name
        inputs = json.loads((folder / "input.json").read_text())
        r_file, py_file = folder / f"screen-{A}.h5mu", folder / f"screen-{B}.h5mu"
        with h5py.File(r_file) as a, h5py.File(py_file) as b:
            for modality in ("rna", "grna"):
                arrays = all(a[f"mod/{modality}/X/{k}"].dtype == b[f"mod/{modality}/X/{k}"].dtype and np.array_equal(a[f"mod/{modality}/X/{k}"][()], b[f"mod/{modality}/X/{k}"][()])
                             for k in ("data", "indices", "indptr"))
                results.add(f"{seed}: {modality} count arrays identical in {A}'s file and {B}'s", bool(arrays))
                for table in ("obs", "var"):
                    verdict = same_tables(a[f"mod/{modality}/{table}"], b[f"mod/{modality}/{table}"])
                    results.add(f"{seed}: {modality} {table} table has the same columns, encodings and values", verdict is True, "" if verdict is True else verdict)
            verdict = same_tables(a["obs"], b["obs"])
            results.add(f"{seed}: shared cell table has the same columns, encodings and values", verdict is True, "" if verdict is True else verdict)
            stamp = lambda f: {k: (f["uns/sceptre"][k].asstr()[()], str(f["uns/sceptre"][k].dtype)) for k in f["uns/sceptre"]}
            results.add(f"{seed}: file metadata identical", stamp(a) == stamp(b), f"{stamp(a)} / {stamp(b)}")
        for language, file in ((A, r_file), (B, py_file)):
            problems = sf.validate_screen(file)
            results.add(f"{seed}: Python validator passes the file {language} wrote", len(problems) == 0, problems.to_string() if len(problems) else "")
            problems = sf.validate_analysis(folder / f"analysis-{language}", file)
            results.add(f"{seed}: Python validator passes the analysis {language} wrote", len(problems) == 0, problems.to_string() if len(problems) else "")

        # the analysis folders, table by table
        r_folder, py_folder = folder / f"analysis-{A}", folder / f"analysis-{B}"
        r_tables = sorted(str(p.relative_to(r_folder)) for p in r_folder.rglob("*.parquet") if "fits" not in p.parts)
        py_tables = sorted(str(p.relative_to(py_folder)) for p in py_folder.rglob("*.parquet") if "fits" not in p.parts)
        results.add(f"{seed}: the two analyses have the same tables", r_tables == py_tables, f"{r_tables} / {py_tables}")
        for name in r_tables:
            x, y = pq.read_table(r_folder / name), pq.read_table(py_folder / name)
            same = x.column_names == y.column_names and [str(t) for t in x.schema.types] == [str(t) for t in y.schema.types] and x.to_pylist() == y.to_pylist()
            results.add(f"{seed}: {name} has the same columns, types and rows", bool(same), f"{x.schema} / {y.schema}" if not same else "")
        fits = lambda f: [(r.response_id, [float(v) for v in r.coefficients], float(r.theta)) for r in sf.read_fits(f)["table"].itertuples()]
        results.add(f"{seed}: fit cache has the same rows, to the last bit", fits(r_folder) == fits(py_folder))
        for stage in sf.STAGES:
            if (r_folder / stage / "record.json").exists():
                strip = lambda rec: {k: (v if k != "provenance" else {**v, "language": ""}) for k, v in rec.items()}
                x, y = json.loads((r_folder / stage / "record.json").read_text()), json.loads((py_folder / stage / "record.json").read_text())
                results.add(f"{seed}: {stage} record equal apart from the writer's language", strip(x) == strip(y), f"{x} / {y}" if strip(x) != strip(y) else "")

        # the first implementation's file read back and compared with the source
        back = md.read_h5mu(r_file)
        for modality, key in (("rna", "rna"), ("grna", "grna")):
            coo = back[modality].X.tocoo()
            read = sorted([int(r), int(c), int(v)] for r, c, v in zip(coo.row, coo.col, coo.data))
            results.add(f"{seed}: {modality} counts read from {A}'s file equal the source", read == sorted(inputs[key]))
        results.add(f"{seed}: identifiers read from {A}'s file equal the source",
                    [list(back["rna"].obs_names), list(back["rna"].var_names), list(back["grna"].var_names)] == [inputs["cells"], inputs["genes"], inputs["grnas"]])
        batch = back.obs["batch"]
        results.add(f"{seed}: factor levels and values read from {A}'s file equal the source",
                    list(batch.cat.categories) == inputs["batch"]["levels"] and list(batch) == inputs["batch"]["values"])
        quality = [None if np.isnan(v) else float(v) for v in back.obs["quality"]]
        results.add(f"{seed}: doubles with missing values read from {A}'s file equal the source", quality == inputs["quality"])
        starts = [None if pd.isna(v) else int(v) for v in back["grna"].var["target_start"]]
        results.add(f"{seed}: gRNA table read from {A}'s file equals the source",
                    list(back["grna"].var["grna_target"]) == inputs["grna_target"] and starts == inputs["target_start"]
                    and [None if v is None or pd.isna(v) else v for v in back["grna"].var["target_chr"]] == inputs["target_chr"])
        saved = sf.read_analysis(r_folder)
        got = saved["tables"]["discovery/results"]
        want = pd.DataFrame(inputs["results"]).sort_values(["grna_target", "response_id"], key=lambda s: s.map(lambda v: v.encode())).reset_index(drop=True)
        equal = all(
            [None if pd.isna(v) else v for v in got[c]] == [None if v is None or (isinstance(v, float) and np.isnan(v)) else v for v in want[c]]
            for c in ("grna_target", "response_id", "pass_qc", "p_value", "fold_change", "se_fold_change"))
        results.add(f"{seed}: results read from {A}'s analysis equal the source, to the last bit", bool(equal))
        results.add(f"{seed}: assignments read from {A}'s analysis equal the source",
                    sorted(map(list, saved["tables"]["assignment/assignments"][["grna_id", "cell_id"]].itertuples(index=False))) == sorted(inputs["assignments"]))
    results.save()
