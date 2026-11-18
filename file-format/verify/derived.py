"""Derived values against frozen expectations. What a reader computes from the
files (which cells each gRNA and each target has, which cells are in use, which
pairs pass QC, the log2 fold change, which results are significant, the control
groups) is compared with what sceptre itself held, stored in
example/<dataset>/expected/. Running this loads neither sceptre nor pysceptre.

    python verify/derived.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "verify")
sys.path.insert(0, "format")
import sceptre_format as sf
from common import Results

results = Results("derived")


def ids_sha(lists: dict, cells) -> str:
    text = "\n".join(f"{name}\t{','.join(cells[np.asarray(idx, dtype=int)])}" for name, idx in lists.items())
    return hashlib.sha256(text.encode()).hexdigest()


for expected_file in sorted(Path("example").glob("*/expected/*/expected.json")):
    dataset, name = expected_file.parts[1], expected_file.parts[3]
    label = f"{dataset}/{name}"
    expected = json.loads(expected_file.read_text())
    screen = sf.read_screen(f"example/{dataset}/screen.h5mu")
    saved = sf.read_analysis(f"example/{dataset}/analyses/{name}")
    tables, records = saved["tables"], saved["stages"]
    cells = np.asarray(screen["grna"].obs_names, dtype=object)
    strategy = records["parameters"]["settings"]["grna_integration_strategy"]
    grna_table = screen["grna"].var.reset_index(names="grna_id")
    grna_table["grna_group"] = sf.grna_groups(grna_table, strategy)
    if "cells_by_grna_sha256" in expected:
        by_grna = sf.cells_by_grna(tables["assignment/assignments"], list(screen["grna"].var_names), pd.Index(cells))
        by_group, by_control = sf.group_cells(by_grna, grna_table)
        results.add(f"{label}: cells of each gRNA", ids_sha(by_grna, cells) == expected["cells_by_grna_sha256"])
        results.add(f"{label}: cells of each target (pooled over its gRNAs)", ids_sha(by_group, cells) == expected["cells_by_group_sha256"])
        results.add(f"{label}: cells of each non-targeting gRNA", ids_sha(by_control, cells) == expected["cells_by_control_grna_sha256"])
        if "grnas_per_cell_sha256" in expected:
            per_cell = ",".join(map(str, sf.grnas_per_cell(by_grna, len(cells))))
            results.add(f"{label}: gRNAs per cell", hashlib.sha256(per_cell.encode()).hexdigest() == expected["grnas_per_cell_sha256"])
    if "cells_in_use_sha256" in expected:
        kept = tables["qc/cells"].set_index("cell_id").loc[cells, "in_use"].to_numpy(dtype=bool)
        results.add(f"{label}: cells in use", sf.sha256_ids(cells[kept]) == expected["cells_in_use_sha256"])
        for pair_type in ("discovery", "positive_control"):
            want = expected.get(f"n_ok_{pair_type}_pairs")
            if isinstance(want, (int, float)):
                results.add(f"{label}: {pair_type} pairs passing QC", sf.n_pairs_passing(tables[f"qc/{pair_type}_pair_qc"], strategy) == want)
    if "control_groups" in expected:
        groups = tables["calibration/negative_control_groups"].groupby("grna_group")["grna_id"].apply(set)
        results.add(f"{label}: the control groups sceptre formed", sorted(groups.index, key=str.encode) == expected["control_groups"])
        non_targeting = set(grna_table.loc[grna_table["grna_target"] == "non-targeting", "grna_id"])
        results.add(f"{label}: every control group's members are non-targeting gRNAs, and are the ones its name lists",
                    all(members <= non_targeting and members == set(group.split("&")) for group, members in groups.items()))
    for stage in ("calibration", "power", "discovery"):
        path = expected_file.parent / f"{stage}.parquet"
        if not path.exists():
            continue
        want = sf.read_table(path)
        got = sf.add_log_2_fold_change(tables[f"{stage}/results"])
        if stage != "power":
            got = sf.add_significant(got, records[stage].get("significance_threshold"))
        key = sf.unit_key(want)
        both = want.merge(got, on=key, suffixes=("_expected", ""))
        results.add(f"{label}: {stage} results cover the same pairs", len(both) == len(want) == len(got))
        for column in ("p_value", "log_2_fold_change"):
            results.add(f"{label}: {stage} {column} identical to sceptre's, to the last bit",
                        np.array_equal(both[column].to_numpy(float), both[f"{column}_expected"].to_numpy(float), equal_nan=True))
        if "significant_expected" in both:
            same = [(None if pd.isna(a) else bool(a)) == (None if pd.isna(b) else bool(b)) for a, b in zip(both["significant"], both["significant_expected"])]
            results.add(f"{label}: {stage} `significant`, from the stored threshold, identical to sceptre's", all(same))
results.save()
