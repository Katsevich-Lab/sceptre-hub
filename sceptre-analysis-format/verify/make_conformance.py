"""Builds conformance/: analysis folders that conform to SPEC.md and analysis
folders that each break one of its rules, with the rules a validator must report
for each (conformance/cases.json). It builds from the committed analyses in
example/, which are fixed: nothing here regenerates them. Each case names the
data file its analysis is checked against, in the data file format's
conformance set. Run from the folder that holds format/, after the data file
format's own conformance set has been built:

    python verify/make_conformance.py

Every invalid folder is a copy of a valid one with one thing changed.
"""

import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import mudata as md
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse

sys.path.insert(0, "format")
import sceptre_analysis as sf

ROOT = Path("conformance")
shutil.rmtree(ROOT, ignore_errors=True)
VALID, INVALID = ROOT / "valid", ROOT / "invalid"
cases = []


def case(kind, path, expect, what, screen=None):
    """`expect` is the rules a validator must report; an analysis conforms when there
    are none. `screen` is the data file the analysis is checked against."""
    entry = {"name": str(Path(path).relative_to(ROOT)), "kind": kind, "path": str(path), "expect": sorted(expect), "conforms": not expect, "what": what}
    if screen is not None:
        entry["screen"] = str(screen)
    cases.append(entry)


# the data files the analyses are checked against
DATA = Path("../perturb-seq-data-format/conformance/valid")
high_r, high_py = DATA / "highmoi_written_by_r.h5mu", DATA / "highmoi_written_by_python.h5mu"
VALID.mkdir(parents=True)

# ---------------------------------------------------------------- valid analyses

for name in ("left_sided", "bonferroni", "singleton", "stopped_after_qc", "with_resamples", "mixture_assignment"):
    shutil.copytree(f"example/highmoi/analyses/{name}", VALID / "analyses" / f"{name}_written_by_r")
    case("analysis", VALID / "analyses" / f"{name}_written_by_r", [], f"The `{name}` example analysis, written by R", high_r)

# the same analysis written again by Python, stage by stage, with the plain reader and writer
saved = sf.read_analysis("example/highmoi/analyses/left_sided")
left_py = VALID / "analyses" / "left_sided_written_by_python"
sf.init_analysis(left_py, saved["screen"], saved["response_modality"])
for stage, record in saved["stages"].items():
    tables = {f"{stage}/{f['path'].removesuffix('.parquet')}": saved["tables"][f"{stage}/{f['path'].removesuffix('.parquet')}"] for f in record["files"]}
    body = {k: v for k, v in record.items() if k not in ("stage", "provenance", "files")}
    sf.write_stage(left_py, stage, body, tables, record["provenance"])
case("analysis", left_py, [], "The `left_sided` analysis, written by Python", high_r)
case("analysis", left_py, [], "The same, checked against the Python-written data file", high_py)
cases[-1]["name"] += " (Python data file)"

parts = VALID / "analyses" / "tables_in_parts"
shutil.copytree("example/highmoi/analyses/left_sided", parts)
table = pq.read_table(parts / "assignment/assignments.parquet")
(parts / "assignment/assignments.parquet").unlink()
(parts / "assignment/assignments.parquet").mkdir()
size = -(-table.num_rows // 3)
for i in range(3):
    part = table.slice(i * size, size)
    part = pa.table({"grna_id": part["grna_id"].dictionary_encode(), "cell_id": part["cell_id"].dictionary_encode(),
                     "score": np.linspace(0, 1, part.num_rows)})
    pq.write_table(part, parts / "assignment/assignments.parquet" / f"part-{i:05d}.parquet")
(parts / "notes.txt").write_text("A file the format does not define.\n")

for name, change, what in (
    ("data_file_at_a_url", lambda v: v["screen"].__setitem__("path", "https://example.org/data/screen.h5mu"), "The hint to the data file is a URL"),
    ("data_file_beside_the_analysis", lambda v: v["screen"].__setitem__("path", "../screen.h5mu"), "The hint to the data file is a relative path"),
):
    folder = VALID / "analyses" / name
    shutil.copytree("example/highmoi/analyses/left_sided", folder)
    value = json.loads((folder / "analysis.json").read_text())
    change(value)
    (folder / "analysis.json").write_text(json.dumps(value, indent=2) + "\n")
    case("analysis", folder, [], what, high_r)

# the fit cache filled by two jobs: two part files, one response in both
two_parts = VALID / "analyses" / "fits_from_two_jobs"
shutil.copytree("example/highmoi/analyses/left_sided", two_parts)
first = next((two_parts / "fits/response").glob("*.parquet"))
fits_table = pq.read_table(first)
fits_made = {sf.FIT_KEY.encode(): pq.ParquetFile(first).metadata.metadata[sf.FIT_KEY.encode()]}
pq.write_table(fits_table.slice(0, 60).replace_schema_metadata(fits_made), two_parts / "fits/response/part-job1.parquet")
pq.write_table(fits_table.slice(59, 41).replace_schema_metadata(fits_made), two_parts / "fits/response/part-job2.parquet")
first.unlink()
case("analysis", two_parts, [], "The fit cache was filled by two jobs, and one response is in both parts", high_r)

# the same formula, spaced differently in the record of how the fits were made
spaced = VALID / "analyses" / "fits_formula_spaced_differently"
shutil.copytree("example/highmoi/analyses/left_sided", spaced)
FIT_KEY = sf.FIT_KEY.encode()


def fit_part(folder):
    return next((Path(folder) / "fits/response").glob("*.parquet"))


def restamp(part, change, to=None):
    """Write a part again with its identity changed by `change` (None removes it)."""
    table = pq.read_table(part)
    # the file's own key-value metadata: pyarrow's table shows only what Arrow's schema
    # record holds, which is not where a part written by R keeps the key
    stored = {k: v for k, v in (pq.ParquetFile(part).metadata.metadata or {}).items() if k != b"ARROW:schema"}
    if change is None:
        stored.pop(FIT_KEY, None)
    else:
        identity = json.loads(stored[FIT_KEY])
        out = change(identity)
        stored[FIT_KEY] = (out if isinstance(out, str) else json.dumps(identity)).encode()
    pq.write_table(table.replace_schema_metadata(stored), to or part)


restamp(fit_part(spaced), lambda v: v.__setitem__("formula", v["formula"].replace(" + ", "+").replace("~", "~ ")))
case("analysis", spaced, [], "A part of the fit cache records the analysis's formula with other spacing; its fits are used", high_r)
cases[-1]["fits_used"] = True

# Parts that are not this analysis's are allowed and unused: nothing depends on
# them, so the validators say nothing and a reader takes no fit from them.
other_method = VALID / "analyses" / "fits_made_for_another_control_group"
shutil.copytree("example/highmoi/analyses/left_sided", other_method)
restamp(fit_part(other_method), lambda v: v.__setitem__("control_group", "nt_cells"))
case("analysis", other_method, [], "The one part of the fit cache was made for a control group other than the analysis's: allowed, and a reader takes no fit from it", high_r)
cases[-1]["fits_used"] = False

stale = VALID / "analyses" / "stale_part_beside_a_current_one"
shutil.copytree("example/highmoi/analyses/left_sided", stale)
current_identity = json.loads(pq.ParquetFile(fit_part(stale)).metadata.metadata[FIT_KEY])
restamp(fit_part(stale), lambda v: v.__setitem__("cells_in_use_sha256", "0" * 64), to=stale / "fits/response/part-left-by-an-earlier-qc.parquet")
case("analysis", stale, [], "Beside the current part is one left by an earlier QC, made for other cells: allowed, and a reader given the current identity reads the current part only", high_r)
cases[-1]["fits_for"] = {"identity": current_identity, "rows": pq.read_table(fit_part(other_method)).num_rows}

# a testing stage rerun: its folder is gone, and the fit cache rightly remains
cleared = VALID / "analyses" / "testing_stage_cleared_fits_kept"
shutil.copytree("example/highmoi/analyses/left_sided", cleared)
sf.clear_stage(cleared, "calibration")
case("analysis", cleared, [], "The calibration stage was cleared for a rerun; the fit cache remains", high_r)

other_language = VALID / "analyses" / "stage_written_in_another_language"
shutil.copytree("example/highmoi/analyses/left_sided", other_language)
record = json.loads((other_language / "power/record.json").read_text())
record["provenance"].update({"package": "Sceptre.jl", "language": "Julia"})
(other_language / "power/record.json").write_text(json.dumps(record, indent=2) + "\n")
case("analysis", other_language, [], "A stage written by an implementation in a third language", high_r)
case("analysis", parts, [], "A table as a folder of parts, with dictionary-encoded identifiers, an unknown column, and an unknown file beside it", high_r)

# ---------------------------------------------------------------- invalid analyses


def broken_analysis(name, expect, what, change, source="example/highmoi/analyses/left_sided", reader_refuses=False):
    path = INVALID / "analyses" / name
    shutil.copytree(source, path)
    change(path)
    case("analysis", path, expect, what, high_r)
    if reader_refuses:
        cases[-1]["reader_refuses"] = True   # a reader must stop, not only the validator


def edit_json(path, change):
    value = json.loads(path.read_text())
    change(value)
    path.write_text(json.dumps(value, indent=2) + "\n")


def rewrite(path, change):
    """Replace a table with `change(table)`, sorted as it was unless `change` says otherwise."""
    # the file's own key-value metadata is kept: a part of the fit cache says there how its fits were made
    kept = {k: v for k, v in (pq.ParquetFile(path).metadata.metadata or {}).items() if k != b"ARROW:schema"}
    pq.write_table(change(pq.read_table(path)).replace_schema_metadata(kept or None), path)


broken_analysis("no_analysis_json", ["A1.1"], "`analysis.json` is missing", lambda p: (p / "analysis.json").unlink())
broken_analysis("no_response_modality", ["A1.1"], "`analysis.json` does not name a response modality",
                lambda p: edit_json(p / "analysis.json", lambda v: v.pop("response_modality")))
broken_analysis("response_modality_grna", ["A1.1"], "`analysis.json` names `grna` as the response modality",
                lambda p: edit_json(p / "analysis.json", lambda v: v.__setitem__("response_modality", "grna")))
broken_analysis("response_modality_not_in_file", ["A1.2"], "The response modality is `protein`, which the data file does not have",
                lambda p: edit_json(p / "analysis.json", lambda v: v.__setitem__("response_modality", "protein")))
# covariates an analysis brings itself (A14.6 to A14.8): here the count summaries of a
# data file that does not store them, and a label that some method computed
BARE = DATA / "highmoi_without_count_summaries.h5mu"


def with_covariates(folder, change=None):
    """The left_sided analysis with its own covariate table: the four count
    summaries its formula names, taken from the data file that has them, and a label."""
    shutil.copytree("example/highmoi/analyses/left_sided", folder)
    with h5py.File(high_py) as f:
        cells = [v.decode() for v in f["mod/rna/obs/_index"][()]]
        table = {"cell_id": cells}
        for modality in ("rna", "grna"):
            for name in ("n_nonzero", "n_umis"):
                table[f"{modality}_{name}"] = f[f"mod/{modality}/obs/{name}"][()].astype(float)
    table["cell_cycle"] = [("G1", "S", "G2M")[i % 3] for i in range(len(cells))]
    table = pa.table(table)
    table = table.take(pa.array(sorted(range(len(cells)), key=lambda i: cells[i].encode())))
    if change is not None:
        table = change(table)
    pq.write_table(table, folder / "parameters/cell_covariates.parquet")
    edit_json(folder / "parameters/record.json", lambda v: v["files"].append({"path": "cell_covariates.parquet", "rows": table.num_rows}))


brought = VALID / "analyses" / "covariates_brought_by_the_analysis"
with_covariates(brought)
case("analysis", brought, [], "The data file stores no count summaries; the analysis brings them, and a computed label, in its own covariate table", BARE)

named_twice = INVALID / "analyses" / "covariate_the_data_file_already_has"
with_covariates(named_twice)
case("analysis", named_twice, ["A14.7"], "The analysis brings covariates under names the data file already provides", high_r)

unknown = INVALID / "analyses" / "formula_names_an_unknown_covariate"
shutil.copytree("example/highmoi/analyses/left_sided", unknown)
case("analysis", unknown, ["A14.8"], "The formula names count summaries that neither the data file nor the analysis provides", BARE)

wrong_type = INVALID / "analyses" / "covariate_of_the_wrong_type"
with_covariates(wrong_type, lambda t: t.append_column("flag", pa.array([True] * t.num_rows)))
case("analysis", wrong_type, ["A14.6"], "A covariate the analysis brings is a boolean, not a number or text", BARE)

too_few = INVALID / "analyses" / "covariates_for_some_cells_only"
with_covariates(too_few, lambda t: t.slice(0, t.num_rows - 10))
case("analysis", too_few, ["A14.6"], "The analysis's covariate table lacks some cells of the data file", BARE)

broken_analysis("analysis_format_from_the_future", ["A1.1"], "The analysis's format version is 9.0",
                lambda p: edit_json(p / "analysis.json", lambda v: v.__setitem__("format", "9.0")))
broken_analysis("absolute_path_to_the_data", ["A1.4"], "The hint to the data file is an absolute local path",
                lambda p: edit_json(p / "analysis.json", lambda v: v["screen"].__setitem__("path", "/home/someone/data/screen.h5mu")))
broken_analysis("another_dataset", ["A1.2"], "The dataset identifier is that of another data file; the parts of the fit cache are then another dataset's, which is allowed and unused",
                lambda p: edit_json(p / "analysis.json", lambda v: v["screen"].__setitem__("id", "0" * 32)))
broken_analysis("record_of_another_stage", ["A2.1"], "The record in `discovery/` says it is the record of the power stage",
                lambda p: edit_json(p / "discovery/record.json", lambda v: v.__setitem__("stage", "power")), reader_refuses=True)
broken_analysis("side_up", ["A2.1"], "The parameters record has `side` set to `up`",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].__setitem__("side", "up")))
broken_analysis("moi_medium", ["A2.1"], "The parameters record has `moi` set to `medium`",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].__setitem__("moi", "medium")))
broken_analysis("no_resampling_budget", ["A2.1"], "The parameters record has no `B1`",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].pop("B1")))
broken_analysis("approximation_unknown", ["A2.1"], "The parameters record has `resampling_approximation` set to a value sceptre does not have",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].__setitem__("resampling_approximation", "normal")))
broken_analysis("assignment_method_unknown", ["A2.1"], "The assignment record names a method sceptre does not have",
                lambda p: edit_json(p / "assignment/record.json", lambda v: v["settings"].__setitem__("method", "nearest")))
broken_analysis("no_pairwise_qc_thresholds", ["A2.1"], "The QC record lacks the two thresholds of pairwise QC",
                lambda p: edit_json(p / "qc/record.json", lambda v: [v["settings"].pop(k) for k in ("n_nonzero_trt_thresh", "n_nonzero_cntrl_thresh")]))
broken_analysis("pair_qc_without_counts", ["A7.1"], "A pair-QC table lacks `n_nonzero_trt` (A10.2)",
                lambda p: rewrite(p / "qc/discovery_pair_qc.parquet", lambda t: t.drop_columns(["n_nonzero_trt"])))
broken_analysis("pair_qc_verdict_against_the_counts", ["A10.4"], "Pairs are marked as passing pairwise QC although a count is below its threshold",
                lambda p: edit_json(p / "qc/record.json", lambda v: v["settings"].__setitem__("n_nonzero_trt_thresh", 50)))
broken_analysis("results_without_counts", ["A7.1"], "The discovery results lack `n_nonzero_cntrl` (A12.1)",
                lambda p: rewrite(p / "discovery/results.parquet", lambda t: t.drop_columns(["n_nonzero_cntrl"])))
broken_analysis("stage_without_its_predecessor", ["A2.2", "A4.1"], "The assignment record is missing while later stages remain",
                lambda p: (p / "assignment/record.json").unlink())
broken_analysis("listed_file_missing", ["A5.1"], "A table that a record lists is not in the folder",
                lambda p: (p / "power/results.parquet").unlink())
broken_analysis("table_truncated", ["A5.2", "A9.1"], "`cells.parquet` has fewer rows than its stage wrote",
                lambda p: rewrite(p / "qc/cells.parquet", lambda t: t.slice(0, 450)))


broken_analysis("table_listed_outside_the_stage_folder", ["A2.1"], "A record lists a table by a path that leaves the stage's folder",
                lambda p: edit_json(p / "qc/record.json", lambda v: v["files"].append({"path": "../parameters/discovery_pairs.parquet", "rows": 1406})),
                reader_refuses=True)
broken_analysis("required_table_not_listed", ["A4.1", "A5.4"], "The QC record does not list `cells.parquet`, which is still in the folder",
                lambda p: edit_json(p / "qc/record.json", lambda v: v.__setitem__("files", [f for f in v["files"] if f["path"] != "cells.parquet"])))
broken_analysis("stale_results", ["A4.1"], "The discovery record is gone and its results table remains",
                lambda p: (p / "discovery/record.json").unlink())
broken_analysis("column_missing", ["A7.1"], "`cells.parquet` has no `in_use` column",
                lambda p: rewrite(p / "qc/cells.parquet", lambda t: t.select(["cell_id"])))
broken_analysis("column_of_the_wrong_type", ["A7.1"], "`pass_qc` is stored as integers",
                lambda p: rewrite(p / "qc/discovery_pair_qc.parquet", lambda t: t.set_column(t.column_names.index("pass_qc"), "pass_qc", t["pass_qc"].cast(pa.int32()))))
broken_analysis("rows_not_sorted", ["A7.2"], "The assignment table is in reverse order",
                lambda p: rewrite(p / "assignment/assignments.parquet", lambda t: t.take(list(range(t.num_rows - 1, -1, -1)))))
broken_analysis("duplicate_row", ["A7.4"], "The assignment table has the same (gRNA, cell) twice",
                lambda p: rewrite(p / "assignment/assignments.parquet", lambda t: t.take([0] + list(range(t.num_rows - 1)))))


def unknown_cell(t):
    frame = t.to_pandas()
    frame.loc[len(frame) - 1, "cell_id"] = "not_a_cell"
    order = sf._byte_order(frame, ["grna_id", "cell_id"])
    return pa.Table.from_pandas(frame.iloc[order].reset_index(drop=True), preserve_index=False)


broken_analysis("cell_not_in_the_data_file", ["A7.5"], "The assignment table names a cell the data file does not have",
                lambda p: rewrite(p / "assignment/assignments.parquet", unknown_cell))
def with_grna_id(t):
    return t.append_column("grna_id", pa.array([f"{g}_grna" for g in t["grna_target"].to_pylist()], type=pa.string()))


broken_analysis("grna_column_under_the_union_strategy", ["A10.3"], "A pair-QC table has a `grna_id` column in an analysis that pools gRNAs by target",
                lambda p: rewrite(p / "qc/discovery_pair_qc.parquet", with_grna_id))
broken_analysis("grna_column_in_bonferroni_discovery", ["A12.6"], "The discovery results of a bonferroni analysis have a `grna_id` column",
                lambda p: rewrite(p / "discovery/results.parquet", with_grna_id), source="example/highmoi/analyses/bonferroni")
broken_analysis("significant_stored", ["A12.3"], "The discovery results store a `significant` column",
                lambda p: rewrite(p / "discovery/results.parquet", lambda t: t.append_column("significant", pa.array([False] * t.num_rows))))
broken_analysis("no_log_2_fold_change", ["A7.1"], "The discovery results have no `log_2_fold_change` (A12.1)",
                lambda p: rewrite(p / "discovery/results.parquet", lambda t: t.drop_columns(["log_2_fold_change"])))
broken_analysis("effect_sizes_disagree", ["A12.2"], "A `log_2_fold_change` is not the base-2 logarithm of the row's `fold_change`",
                lambda p: rewrite(p / "discovery/results.parquet", lambda t: t.set_column(
                    t.schema.get_field_index("log_2_fold_change"), "log_2_fold_change", pa.array([None if v is None else v + 1.0 for v in t["log_2_fold_change"].to_pylist()], type=pa.float64()))))
broken_analysis("fit_part_without_identity", ["A13.1"], "A part of the fit cache does not say how its fits were made, as after a tool rewrote it; a reader takes no fit from it",
                lambda p: restamp(fit_part(p), None))
cases[-1]["fits_used"] = False
broken_analysis("fit_identity_incomplete", ["A13.1"], "The identity inside a part of the fit cache lacks `cells_in_use_sha256`",
                lambda p: restamp(fit_part(p), lambda v: v.pop("cells_in_use_sha256") and None))
broken_analysis("fit_identity_not_json", ["A13.1"], "The identity inside a part of the fit cache is not JSON",
                lambda p: restamp(fit_part(p), lambda v: "poisson, all cells"))
cases[-1]["fits_used"] = False
broken_analysis("identity_file_in_the_fit_cache", ["A13.1"], "The fit cache holds a JSON file beside its parts, which a folder reader would trip on",
                lambda p: (p / "fits/response/identity.json").write_text("{}\n"))
broken_analysis("fits_before_qc", ["A13.3"], "An analysis that stopped after assignment holds a fit cache",
                lambda p: [shutil.rmtree(p / s) for s in ("qc", "calibration", "power", "discovery")])
broken_analysis("fit_with_too_few_coefficients", ["A13.2"], "A saved fit has fewer coefficients than the design has terms",
                lambda p: rewrite(next((p / "fits/response").glob("*.parquet")), lambda t: t.set_column(
                    t.column_names.index("coefficients"), "coefficients",
                    pa.array([c[:-1] for c in t["coefficients"].to_pylist()], type=pa.list_(pa.float64())))))
broken_analysis("formula_with_an_interaction", ["A14.1"], "The formula has an interaction, outside the additive subset; the fit cache was made for another formula, which is allowed and unused",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].__setitem__("formula", "~log(rna_n_umis) * batch")))

(ROOT / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
conforming = sum(1 for c in cases if c["conforms"])
print(f"{len(cases)} cases: {conforming} conform, {len(cases) - conforming} do not; "
      f"rules covered: {len({r for c in cases for r in c['expect']})}")
