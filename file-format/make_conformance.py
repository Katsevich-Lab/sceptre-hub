"""Builds conformance/: files that conform to SPEC.md and files that each break
one of its rules, with the rules a validator must report for each
(conformance/cases.json). Run from this folder after make_example.R:

    python make_conformance.py

The valid files are written by R (copied from example/) and by Python (written
here with sceptre_format.py). Every invalid file is a copy of a valid one with
one thing changed.
"""

import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, "format")
import sceptre_format as sf

ROOT = Path("conformance")
shutil.rmtree(ROOT, ignore_errors=True)
VALID, INVALID = ROOT / "valid", ROOT / "invalid"
cases = []


def case(kind, path, expect, what, screen=None, reads_as=None):
    """`level` is what a validator must conclude: "conforming"; "readable" (not
    conforming, and a reader must still accept the file, giving the same data as
    the file `reads_as`); or "not readable"."""
    level = "conforming" if not expect else ("readable" if reads_as is not None else "not readable")
    entry = {"name": str(Path(path).relative_to(ROOT)), "kind": kind, "path": str(path), "expect": sorted(expect), "level": level, "what": what}
    if screen is not None:
        entry["screen"] = str(screen)
    if reads_as is not None:
        entry["reads_as"] = str(reads_as)
    cases.append(entry)


# ---------------------------------------------------------------- valid files

types_r = VALID / "types_written_by_r.h5mu"
types_py = VALID / "types_written_by_python.h5mu"
high_r = VALID / "highmoi_written_by_r.h5mu"
high_py = VALID / "highmoi_written_by_python.h5mu"
VALID.mkdir(parents=True)
shutil.copy("example/types/screen.h5mu", types_r)
shutil.copy("example/highmoi/screen.h5mu", high_r)
sf.write_screen(sf.read_screen(types_r), types_py)
sf.write_screen(sf.read_screen(high_r), high_py)
case("screen", types_r, [], "Every allowed kind of column, written by R")
case("screen", types_py, [], "The same data, written by Python")
case("screen", high_r, [], "The high-MOI example, written by R")
case("screen", high_py, [], "The high-MOI example, written by Python")

gzip = VALID / "types_gzip_counts.h5mu"
shutil.copy(types_r, gzip)
with h5py.File(gzip, "r+") as f:
    for name in ("data", "indices", "indptr"):
        values, attributes = f[f"mod/rna/X/{name}"][()], dict(f[f"mod/rna/X/{name}"].attrs)
        del f[f"mod/rna/X/{name}"]
        f["mod/rna/X"].create_dataset(name, data=values, compression="gzip")
        f[f"mod/rna/X/{name}"].attrs.update(attributes)
case("screen", gzip, [], "Counts compressed with gzip, which S11 allows")

extra = VALID / "types_extra_modality_and_columns.h5mu"
import anndata as ad
import mudata as md
from scipy import sparse

base = sf.read_screen(types_py)
hashing = ad.AnnData(
    X=sparse.csc_matrix(np.array([[5, 0], [0, 7], [3, 0], [0, 2], [4, 0], [0, 6]])),
    obs=pd.DataFrame(index=base["grna"].obs_names), var=pd.DataFrame(index=["hto_1", "hto_2"]),
)
grna_var = base["grna"].var.assign(spacer=["ACGT", "TTGA"])
with_extra = md.MuData({"rna": base["rna"], "grna": ad.AnnData(X=base["grna"].X, obs=base["grna"].obs, var=grna_var), "hashing": hashing})
for name in base.obs.columns:
    with_extra.obs[name] = base.obs[name].array
with_extra.uns["sceptre"] = dict(base.uns["sceptre"])
sf.write_screen(with_extra, extra)
case("screen", extra, [], "A modality no analysis uses, and a gRNA column the format does not define")
cases[-1]["copy_keeps_modalities"] = "rna"   # read with this response modality and written again, every modality is still there

# several gRNAs on one vector: each target's gRNAs share a vector, and each
# non-targeting gRNA is alone on a vector named after it
def vectors_of(f):
    ids = [v.decode() for v in f["mod/grna/var/_index"][()]]
    node = f["mod/grna/var/grna_target"]   # plain text, or a categorical
    if isinstance(node, h5py.Group):
        categories = [v.decode() for v in node["categories"][()]]
        targets = [categories[c] for c in node["codes"][()]]
    else:
        targets = [v.decode() for v in node[()]]
    return ids, targets, [i if t == "non-targeting" else f"{t}_vector" for i, t in zip(ids, targets)]


def add_vectors(f, vectors):
    table = f["mod/grna/var"]
    column = table.create_dataset("vector_id", data=np.array(vectors, dtype=object), dtype=h5py.string_dtype())
    column.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
    # readers take a table's columns from this attribute
    order = [c.decode() if isinstance(c, bytes) else c for c in table.attrs["column-order"]]
    table.attrs["column-order"] = np.array(order + ["vector_id"], dtype=h5py.string_dtype())


with_vectors = VALID / "highmoi_grnas_on_shared_vectors.h5mu"
shutil.copy(high_py, with_vectors)
with h5py.File(with_vectors, "r+") as f:
    add_vectors(f, vectors_of(f)[2])
case("screen", with_vectors, [], "Each target's two gRNAs are on one vector (`vector_id`); a non-targeting gRNA's vector is named after it")

for name in ("left_sided", "bonferroni", "singleton", "stopped_after_qc", "with_resamples"):
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
    ("no_hint_to_the_data_file", lambda v: v["screen"].pop("path"), "The analysis gives no hint of where its data file is"),
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
pq.write_table(fits_table.slice(0, 60), two_parts / "fits/response/part-job1.parquet")
pq.write_table(fits_table.slice(59, 41), two_parts / "fits/response/part-job2.parquet")
first.unlink()
case("analysis", two_parts, [], "The fit cache was filled by two jobs, and one response is in both parts", high_r)

# the same formula, spaced differently in the record of how the fits were made
spaced = VALID / "analyses" / "fits_formula_spaced_differently"
shutil.copytree("example/highmoi/analyses/left_sided", spaced)
identity = json.loads((spaced / "fits/response/identity.json").read_text())
identity["formula"] = identity["formula"].replace(" + ", "+").replace("~", "~ ")
(spaced / "fits/response/identity.json").write_text(json.dumps(identity, indent=2) + "\n")
case("analysis", spaced, [], "The fit cache records the analysis's formula with other spacing", high_r)

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

# ---------------------------------------------------------------- invalid data files


def broken_screen(name, expect, what, change, source=types_py, readable=False, reader_refuses=False):
    path = INVALID / "screens" / f"{name}.h5mu"
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(source, path)
    with h5py.File(path, "r+") as f:
        change(f)
    case("screen", path, expect, what, reads_as=source if readable else None)
    if reader_refuses:
        cases[-1]["reader_refuses"] = True   # a reader must stop, not only the validator


def replace(f, path, values, dtype=None):
    attributes = dict(f[path].attrs)
    del f[path]
    f.create_dataset(path, data=values, dtype=dtype)
    f[path].attrs.update(attributes)


def add_column(f, table, name, values, encoding="array", dtype=None):
    column = f[table].create_dataset(name, data=values, dtype=dtype)
    column.attrs.update({"encoding-type": encoding, "encoding-version": "0.2.0"})


text = h5py.string_dtype()
path = INVALID / "screens" / "not_hdf5.h5mu"
path.parent.mkdir(parents=True)
path.write_text("This is not an HDF5 file.\n")
case("screen", path, ["S1.1"], "A text file with the right extension")

broken_screen("no_grna_modality", ["S2.1"], "The gRNA modality is named `guide`", lambda f: f.move("mod/grna", "mod/guide"))
broken_screen("bad_modality_name", ["S2.2"], "A modality named `1rna`", lambda f: f.move("mod/rna", "mod/1rna"))


def reorder_cells(f):
    ids = f["mod/grna/obs/_index"][()]
    ids[[0, 1]] = ids[[1, 0]]
    replace(f, "mod/grna/obs/_index", ids, text)


broken_screen("cells_in_another_order", ["S2.3"], "The gRNA modality lists two cells in the other order", reorder_cells)
broken_screen("compressed_by_cell", ["S3.1"], "The RNA counts are marked as compressed by cell (`csr_matrix`), as in an IGVF file",
              lambda f: f["mod/rna/X"].attrs.__setitem__("encoding-type", "csr_matrix"))


def bad_pointers(f):
    pointers = f["mod/rna/X/indptr"][()]
    pointers[-1] += 1
    replace(f, "mod/rna/X/indptr", pointers)


broken_screen("pointers_past_the_end", ["S3.2"], "The last column pointer is larger than the number of stored values", bad_pointers)


def unsorted_cells(f):
    positions = f["mod/rna/X/indices"][()]
    start, stop = f["mod/rna/X/indptr"][0:2]
    positions[[start, start + 1]] = positions[[start + 1, start]]
    replace(f, "mod/rna/X/indices", positions)


broken_screen("cells_unsorted_within_a_feature", ["S3.3"], "Two cell positions of the first gene are swapped", unsorted_cells)


def stored_zero(f):
    counts = f["mod/rna/X/data"][()]
    counts[0] = 0
    replace(f, "mod/rna/X/data", counts)


broken_screen("stored_zero", ["S3.4"], "A zero is stored as if it were a count", stored_zero)
broken_screen("counts_as_doubles", ["S4.1"], "Counts stored as doubles, as R's MuData writes them by default",
              lambda f: replace(f, "mod/rna/X/data", f["mod/rna/X/data"][()].astype(np.float64)), readable=True)
broken_screen("counts_as_uint16", ["S4.1"], "Counts stored as 16-bit unsigned integers, as an IGVF file has them",
              lambda f: replace(f, "mod/rna/X/data", f["mod/rna/X/data"][()].astype(np.uint16)), readable=True)
# A small file's index type is 32-bit (S4). Wider arrays, in one or both, are readable.
broken_screen("cell_positions_64_bit", ["S4.2"], "Cell positions stored as 64-bit integers beside 32-bit column pointers",
              lambda f: replace(f, "mod/rna/X/indices", f["mod/rna/X/indices"][()].astype(np.int64)), readable=True)
broken_screen("pointers_64_bit", ["S4.3"], "Column pointers stored as 64-bit integers beside 32-bit cell positions: mixed widths, which scipy rejects in some operations",
              lambda f: replace(f, "mod/rna/X/indptr", f["mod/rna/X/indptr"][()].astype(np.int64)), readable=True)


def both_wide(f):
    for name in ("indices", "indptr"):
        replace(f, f"mod/rna/X/{name}", f[f"mod/rna/X/{name}"][()].astype(np.int64))


broken_screen("both_index_arrays_64_bit", ["S4.2", "S4.3"], "Both index arrays stored as 64-bit integers in a file small enough for 32-bit. "
              "The same file is what a file past the limit must look like: with the limit lowered for the test, the validators pass it",
              both_wide, readable=True)
cases[-1]["conforming_with_index_limit"] = 5   # the rna modality has 6 cells
broken_screen("shared_column_named_like_a_modality", ["S5.1"], "A shared cell column named `rna_quality`",
              lambda f: add_column(f, "obs", "rna_quality", np.arange(6, dtype=np.float64)))
broken_screen("covariate_name_with_a_dash", ["S5.4"], "A cell column named `donor-age`, which no formula can use",
              lambda f: add_column(f, "obs", "donor-age", np.arange(6, dtype=np.float64)))
broken_screen("no_n_nonzero", ["S6.1"], "The gRNA modality's cell table has no `n_nonzero`", lambda f: f["mod/grna/obs"].__delitem__("n_nonzero"))


def n_umis_as_labels(f):
    del f["mod/rna/obs/n_umis"]
    column = f["mod/rna/obs"].create_group("n_umis")
    column.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
    add_column(f, "mod/rna/obs/n_umis", "codes", np.array([0, 1, 0, 1, 0, 1], dtype=np.int8))
    add_column(f, "mod/rna/obs/n_umis", "categories", np.array(["high", "low"], dtype=object), "string-array", text)


broken_screen("n_umis_as_labels", ["S6.2"], "`n_umis` stored as a categorical instead of a number", n_umis_as_labels)
broken_screen("boolean_column", ["S7.1"], "A true/false cell column stored as booleans, as Python writes it by default",
              lambda f: add_column(f, "mod/rna/obs", "treated", np.array([True, False] * 3)))
broken_screen("integer_column_64_bit", ["S7.1"], "A 64-bit integer cell column, as pandas writes integers by default",
              lambda f: add_column(f, "mod/rna/obs", "lane_number", np.arange(6, dtype=np.int64)))
broken_screen("text_covariate", ["S7.1"], "A cell covariate stored as plain text instead of a categorical",
              lambda f: add_column(f, "obs", "donor", np.array(["a", "b"] * 3, dtype=object), "string-array", text))
broken_screen("nullable_string_column", ["S7.2"], "A gRNA column in anndata's nullable string encoding, which R cannot open",
              lambda f: f["mod/grna/var"].create_group("note").attrs.update({"encoding-type": "nullable-string-array", "encoding-version": "0.1.0"}))


def bad_codes(f):
    codes = f["obs/batch/codes"][()]
    codes[0] = 7
    replace(f, "obs/batch/codes", codes)


broken_screen("categorical_code_out_of_range", ["S7.3"], "A categorical with a code beyond its categories", bad_codes)
broken_screen("integer_column_with_missing_value", ["S7.4"], "An integer cell column holding R's missing-integer value",
              lambda f: add_column(f, "mod/rna/obs", "plate", np.array([1, 2, -2147483648, 4, 5, 6], dtype=np.int32)))
broken_screen("no_grna_target", ["S8.1"], "The gRNA table has no `grna_target`", lambda f: f["mod/grna/var"].__delitem__("grna_target"))
broken_screen("grna_group_stored", ["S8.5"], "The gRNA table stores `grna_group`, which depends on the analysis",
              lambda f: add_column(f, "mod/grna/var", "grna_group", np.array(["site_1", "non-targeting"], dtype=object), "string-array", text))
def vectors_changed(change):
    def apply(f):
        ids, targets, vectors = vectors_of(f)
        change(ids, targets, vectors)
        add_vectors(f, vectors)
    return apply


def two_targets_on_a_vector(ids, targets, vectors):
    other = next(i for i, t in enumerate(targets) if t not in (targets[0], "non-targeting"))
    vectors[other] = vectors[0]


def vector_named_like_another_grna(ids, targets, vectors):
    other = next(i for i, t in enumerate(targets) if t not in (targets[0], "non-targeting"))
    for i, t in enumerate(targets):
        if t == targets[0]:
            vectors[i] = ids[other]


broken_screen("vector_with_two_targets", ["S8.7"], "Two gRNAs with different targets share a `vector_id`: a combinatorial screen",
              vectors_changed(two_targets_on_a_vector), source=high_py, reader_refuses=True)
broken_screen("vector_named_like_another_grna", ["S8.8"], "A vector is named with the identifier of a gRNA on another vector",
              vectors_changed(vector_named_like_another_grna), source=high_py, reader_refuses=True)
broken_screen("vector_id_as_a_number", ["S8.6"], "`vector_id` stored as numbers",
              lambda f: add_column(f, "mod/grna/var", "vector_id", np.arange(f["mod/grna/var/_index"].shape[0], dtype=np.float64)), source=high_py)
def set_positions(start=None, end=None, chromosome=None):
    """Change the target position of the first gRNA of the types file (chr1, 100 to 149)."""
    def apply(f):
        table = f["mod/grna/var"]
        for name, value in (("target_start", start), ("target_end", end)):
            if value is not None:
                values = table[name][()]
                values[0] = value
                replace(f, f"mod/grna/var/{name}", values)
        if chromosome == "missing":
            codes = table["target_chr/codes"][()]
            codes[0] = -1
            replace(f, "mod/grna/var/target_chr/codes", codes)
    return apply


broken_screen("target_start_zero", ["S8.10"], "A target starts at position 0: a 0-based coordinate, as in a BED file", set_positions(start=0))
broken_screen("target_end_before_start", ["S8.10"], "A target ends before it starts", set_positions(end=50))
broken_screen("target_start_without_chromosome", ["S8.11"], "A gRNA has a target start and end and no chromosome", set_positions(chromosome="missing"))
broken_screen("target_start_as_text", ["S8.2"], "`target_start` stored as text",
              lambda f: (f["mod/grna/var"].__delitem__("target_start"),
                         add_column(f, "mod/grna/var", "target_start", np.array(["100", "200"], dtype=object), "string-array", text)))


def duplicate_cells(f):
    for table in ("obs", "mod/rna/obs", "mod/grna/obs"):
        ids = f[f"{table}/_index"][()]
        ids[1] = ids[0]
        replace(f, f"{table}/_index", ids, text)


broken_screen("duplicate_cell_identifiers", ["S9.1"], "Two cells have the same identifier", duplicate_cells)


def duplicate_features(f):
    ids = f["mod/rna/var/_index"][()]
    ids[0] = f["mod/grna/var/_index"][0]
    replace(f, "mod/rna/var/_index", ids, text)


broken_screen("gene_named_like_a_grna", ["S9.2"], "A gene and a gRNA have the same identifier", duplicate_features)


def line_feed(f):
    for table in ("obs", "mod/rna/obs", "mod/grna/obs"):
        ids = f[f"{table}/_index"][()]
        ids[0] = "cell\n1"
        replace(f, f"{table}/_index", ids, text)


broken_screen("identifier_with_a_line_feed", ["S9.4"], "A cell identifier contains a line feed", line_feed)
broken_screen("no_sceptre_metadata", ["S10.1"], "The file has no `uns/sceptre`", lambda f: f["uns"].__delitem__("sceptre"))


def set_scalar(f, name, value):
    del f[f"uns/sceptre/{name}"]
    f["uns/sceptre"].create_dataset(name, data=value, dtype=text if isinstance(value, str) else None)
    encoding = "string" if isinstance(value, str) else "numeric-scalar"
    f[f"uns/sceptre/{name}"].attrs.update({"encoding-type": encoding, "encoding-version": "0.2.0"})


broken_screen("dataset_name_as_a_number", ["S10.5"], "The optional dataset name is a number",
              lambda f: f["uns/sceptre"].create_dataset("name", data=np.int32(7)) if "name" not in f["uns/sceptre"] else set_scalar(f, "name", np.int32(7)))
broken_screen("unknown_perturbation_type", ["S10.6"], "The optional perturbation type is not one of the listed kinds",
              lambda f: set_scalar(f, "perturbation_type", "CRISPRx"))
broken_screen("moi_medium", ["S10.2"], "`moi` is `medium`", lambda f: set_scalar(f, "moi", "medium"))
broken_screen("format_from_the_future", ["S10.3"], "The data file's format version is 9.0", lambda f: set_scalar(f, "format", "9.0"))
broken_screen("format_as_a_number", ["S10.1"], "The data file's format is the integer 0, not a version", lambda f: set_scalar(f, "format", np.int32(0)))

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
    pq.write_table(change(pq.read_table(path)), path)


broken_analysis("no_analysis_json", ["A1.1"], "`analysis.json` is missing", lambda p: (p / "analysis.json").unlink())
broken_analysis("no_response_modality", ["A1.1"], "`analysis.json` does not name a response modality",
                lambda p: edit_json(p / "analysis.json", lambda v: v.pop("response_modality")))
broken_analysis("response_modality_grna", ["A1.1"], "`analysis.json` names `grna` as the response modality",
                lambda p: edit_json(p / "analysis.json", lambda v: v.__setitem__("response_modality", "grna")))
broken_analysis("response_modality_not_in_file", ["A1.2"], "The response modality is `protein`, which the data file does not have",
                lambda p: edit_json(p / "analysis.json", lambda v: v.__setitem__("response_modality", "protein")))
broken_analysis("analysis_format_from_the_future", ["A1.1"], "The analysis's format version is 9.0",
                lambda p: edit_json(p / "analysis.json", lambda v: v.__setitem__("format", "9.0")))
broken_analysis("absolute_path_to_the_data", ["A1.4"], "The hint to the data file is an absolute local path",
                lambda p: edit_json(p / "analysis.json", lambda v: v["screen"].__setitem__("path", "/home/someone/data/screen.h5mu")))
broken_analysis("another_dataset", ["A1.2", "T6.3"], "The dataset identifier is that of another data file; the fit cache then belongs to another dataset too",
                lambda p: edit_json(p / "analysis.json", lambda v: v["screen"].__setitem__("id", "0" * 32)))
broken_analysis("record_of_another_stage", ["A2.1"], "The record in `discovery/` says it is the record of the power stage",
                lambda p: edit_json(p / "discovery/record.json", lambda v: v.__setitem__("stage", "power")), reader_refuses=True)
broken_analysis("side_up", ["A2.1"], "The parameters record has `side` set to `up`",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].__setitem__("side", "up")))
broken_analysis("no_seed", ["A2.1"], "The parameters record has no `seed`",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].pop("seed")))
broken_analysis("stage_without_its_predecessor", ["A2.2", "A4.1"], "The assignment record is missing while later stages remain",
                lambda p: (p / "assignment/record.json").unlink())
broken_analysis("listed_file_missing", ["A5.1"], "A table that a record lists is not in the folder",
                lambda p: (p / "power/results.parquet").unlink())
broken_analysis("table_truncated", ["A5.2", "T2.1"], "`cells.parquet` has fewer rows than its stage wrote",
                lambda p: rewrite(p / "qc/cells.parquet", lambda t: t.slice(0, 450)))


broken_analysis("table_listed_outside_the_stage_folder", ["A2.1"], "A record lists a table by a path that leaves the stage's folder",
                lambda p: edit_json(p / "qc/record.json", lambda v: v["files"].append({"path": "../parameters/discovery_pairs.parquet", "rows": 1406})),
                reader_refuses=True)
broken_analysis("required_table_not_listed", ["A4.1", "A5.4"], "The QC record does not list `cells.parquet`, which is still in the folder",
                lambda p: edit_json(p / "qc/record.json", lambda v: v.__setitem__("files", [f for f in v["files"] if f["path"] != "cells.parquet"])))
broken_analysis("stale_results", ["A4.1"], "The discovery record is gone and its results table remains",
                lambda p: (p / "discovery/record.json").unlink())
broken_analysis("column_missing", ["T0.1"], "`cells.parquet` has no `in_use` column",
                lambda p: rewrite(p / "qc/cells.parquet", lambda t: t.select(["cell_id"])))
broken_analysis("column_of_the_wrong_type", ["T0.1"], "`pass_qc` is stored as integers",
                lambda p: rewrite(p / "qc/discovery_pair_qc.parquet", lambda t: t.set_column(t.column_names.index("pass_qc"), "pass_qc", t["pass_qc"].cast(pa.int32()))))
broken_analysis("rows_not_sorted", ["T0.2"], "The assignment table is in reverse order",
                lambda p: rewrite(p / "assignment/assignments.parquet", lambda t: t.take(list(range(t.num_rows - 1, -1, -1)))))
broken_analysis("duplicate_row", ["T0.4"], "The assignment table has the same (gRNA, cell) twice",
                lambda p: rewrite(p / "assignment/assignments.parquet", lambda t: t.take([0] + list(range(t.num_rows - 1)))))


def unknown_cell(t):
    frame = t.to_pandas()
    frame.loc[len(frame) - 1, "cell_id"] = "not_a_cell"
    order = sf._byte_order(frame, ["grna_id", "cell_id"])
    return pa.Table.from_pandas(frame.iloc[order].reset_index(drop=True), preserve_index=False)


broken_analysis("cell_not_in_the_data_file", ["T0.5"], "The assignment table names a cell the data file does not have",
                lambda p: rewrite(p / "assignment/assignments.parquet", unknown_cell))
def with_grna_id(t):
    return t.append_column("grna_id", pa.array([f"{g}_grna" for g in t["grna_target"].to_pylist()], type=pa.string()))


broken_analysis("grna_column_under_the_union_strategy", ["T3.3"], "A pair-QC table has a `grna_id` column in an analysis that pools gRNAs by target",
                lambda p: rewrite(p / "qc/discovery_pair_qc.parquet", with_grna_id))
broken_analysis("grna_column_in_bonferroni_discovery", ["T5.6"], "The discovery results of a bonferroni analysis have a `grna_id` column",
                lambda p: rewrite(p / "discovery/results.parquet", with_grna_id), source="example/highmoi/analyses/bonferroni")
broken_analysis("significant_stored", ["T5.3"], "The discovery results store a `significant` column",
                lambda p: rewrite(p / "discovery/results.parquet", lambda t: t.append_column("significant", pa.array([False] * t.num_rows))))
broken_analysis("no_effect_size", ["T5.1"], "The discovery results have neither `fold_change` nor `log_2_fold_change`",
                lambda p: rewrite(p / "discovery/results.parquet", lambda t: t.drop_columns(["fold_change", "se_fold_change"])))
broken_analysis("fits_without_their_record", ["T6.1"], "The fit cache has fits and no record of how they were made",
                lambda p: (p / "fits/response/identity.json").unlink())
broken_analysis("fits_made_with_another_method", ["T6.3"], "The fit cache says its fits were made with a method other than the analysis's",
                lambda p: edit_json(p / "fits/response/identity.json", lambda v: v.__setitem__("response_fit_method", "sceptre")))
broken_analysis("fits_before_qc", ["T6.3"], "An analysis that stopped after assignment holds a fit cache",
                lambda p: [shutil.rmtree(p / s) for s in ("qc", "calibration", "power", "discovery")])
broken_analysis("fit_with_too_few_coefficients", ["T6.2"], "A saved fit has fewer coefficients than the design has terms",
                lambda p: rewrite(next((p / "fits/response").glob("*.parquet")), lambda t: t.set_column(
                    t.column_names.index("coefficients"), "coefficients",
                    pa.array([c[:-1] for c in t["coefficients"].to_pylist()], type=pa.list_(pa.float64())))))
broken_analysis("formula_with_an_interaction", ["T6.3", "T7.1"], "The formula has an interaction, outside the additive subset; the fit cache was made for another formula",
                lambda p: edit_json(p / "parameters/record.json", lambda v: v["settings"].__setitem__("formula", "~log(rna_n_umis) * batch")))

(ROOT / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
levels = {level: sum(1 for c in cases if c["level"] == level) for level in ("conforming", "readable", "not readable")}
print(f"{len(cases)} cases: {levels}; rules covered: {len({r for c in cases for r in c['expect']})}")
