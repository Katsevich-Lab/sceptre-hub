"""Builds conformance/: data files that conform to SPEC.md and data files that
each break one of its rules, with the rules a validator must report for each
(conformance/cases.json). It builds from the committed files in example/, which
are fixed: nothing here regenerates them. Run from the folder that holds format/:

    python verify/make_conformance.py

The valid files are written by R (copied from example/) and by Python (written
here with perturb_seq_data.py). Every invalid file is a copy of a valid one with
one thing changed.
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
from scipy import sparse

sys.path.insert(0, "format")
import perturb_seq_data as sf

ROOT = Path("conformance")
shutil.rmtree(ROOT, ignore_errors=True)
VALID, INVALID = ROOT / "valid", ROOT / "invalid"
cases = []


def case(kind, path, expect, what, reads_as=None):
    """`expect` is the rules a validator must report; a file conforms when there are
    none. `reads_as` names a conforming file that a reader must get the same counts
    from, although this one does not conform (S4.4)."""
    entry = {"name": str(Path(path).relative_to(ROOT)), "kind": kind, "path": str(path), "expect": sorted(expect), "conforms": not expect, "what": what}
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
sf.write_screen(md.read_h5mu(types_r), types_py)
sf.write_screen(md.read_h5mu(high_r), high_py)
case("screen", types_r, [], "Every allowed kind of column, written by R")
case("screen", types_py, [], "The same data, written by Python")
case("screen", high_r, [], "The high-MOI example, written by R")
case("screen", high_py, [], "The high-MOI example, written by Python")

other_tool = VALID / "types_with_another_tools_entry.h5mu"
shutil.copy(types_py, other_tool)
with h5py.File(other_tool, "r+") as f:
    group = f["uns"].create_group("neighbors")
    group.attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
    group.create_dataset("n_neighbors", data=np.int64(15)).attrs.update({"encoding-type": "numeric-scalar", "encoding-version": "0.2.0"})
case("screen", other_tool, [], "`uns` also holds an entry of another tool, which the format ignores (S10.2)")

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
base = md.read_h5mu(types_py)
hashing = ad.AnnData(
    X=sparse.csc_matrix(np.array([[5, 0], [0, 7], [3, 0], [0, 2], [4, 0], [0, 6]])),
    obs=pd.DataFrame(index=base["grna"].obs_names), var=pd.DataFrame(index=["hto_1", "hto_2"]),
)
grna_var = base["grna"].var.assign(spacer=["ACGT", "TTGA"])
with_extra = md.MuData({"rna": base["rna"], "grna": ad.AnnData(X=base["grna"].X, obs=base["grna"].obs, var=grna_var), "hashing": hashing})
for name in base.obs.columns:
    with_extra.obs[name] = base.obs[name].array
with_extra.uns.update({key: base.uns[key] for key in sf.METADATA if key in base.uns})
sf.write_screen(with_extra, extra)
case("screen", extra, [], "A third modality, and a gRNA column the format does not define")
cases[-1]["copy_keeps_modalities"] = "rna"   # read with this modality as R's main experiment and written again, every modality is still there (S2.5)

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


def add_text_column(f, name, values):
    table = f["mod/grna/var"]
    column = table.create_dataset(name, data=np.array(values, dtype=object), dtype=h5py.string_dtype())
    column.attrs.update({"encoding-type": "string-array", "encoding-version": "0.2.0"})
    order = [c.decode() if isinstance(c, bytes) else c for c in table.attrs["column-order"]]
    table.attrs["column-order"] = np.array(order + [name], dtype=h5py.string_dtype())


def types_of(f):
    """The kind of each gRNA of the file: its controls are of two kinds."""
    ids, targets, _ = vectors_of(f)
    controls = iter(["non-targeting", "safe-targeting"] * len(ids))
    return ["targeting" if t != "non-targeting" else next(controls) for t in targets]


def drop_column(f, table, name):
    del f[f"{table}/{name}"]
    order = [c.decode() if isinstance(c, bytes) else c for c in f[table].attrs["column-order"]]
    f[table].attrs["column-order"] = np.array([c for c in order if c != name], dtype=h5py.string_dtype())


without_summaries = VALID / "highmoi_without_count_summaries.h5mu"
shutil.copy(high_py, without_summaries)
with h5py.File(without_summaries, "r+") as f:
    for modality in f["mod"]:
        for name in ("n_umis", "n_nonzero"):
            drop_column(f, f"mod/{modality}/obs", name)
case("screen", without_summaries, [], "No cell table has `n_umis` or `n_nonzero`, which are optional (S6.1)")

with_types = VALID / "highmoi_with_grna_types.h5mu"
shutil.copy(high_py, with_types)
with h5py.File(with_types, "r+") as f:
    add_text_column(f, "grna_type", types_of(f))
case("screen", with_types, [], "The gRNA table says what kind of gRNA each is (`grna_type`); the controls are of two kinds")

with_vectors = VALID / "highmoi_grnas_on_shared_vectors.h5mu"
shutil.copy(high_py, with_vectors)
with h5py.File(with_vectors, "r+") as f:
    add_vectors(f, vectors_of(f)[2])
case("screen", with_vectors, [], "Each target's two gRNAs are on one vector (`vector_id`); a non-targeting gRNA's vector is named after it")

# ---------------------------------------------------------------- invalid data files


def broken_screen(name, expect, what, change, source=types_py, still_read=False, reader_refuses=False):
    path = INVALID / "screens" / f"{name}.h5mu"
    path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(source, path)
    with h5py.File(path, "r+") as f:
        change(f)
    case("screen", path, expect, what, reads_as=source if still_read else None)
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
              lambda f: replace(f, "mod/rna/X/data", f["mod/rna/X/data"][()].astype(np.float64)), still_read=True)
broken_screen("counts_as_uint16", ["S4.1"], "Counts stored as 16-bit unsigned integers, as an IGVF file has them",
              lambda f: replace(f, "mod/rna/X/data", f["mod/rna/X/data"][()].astype(np.uint16)), still_read=True)
# A small file's index type is 32-bit (S4). Wider arrays, in one or both, do not conform, and a reader still reads them (S4.4).
broken_screen("cell_positions_64_bit", ["S4.2"], "Cell positions stored as 64-bit integers beside 32-bit column pointers",
              lambda f: replace(f, "mod/rna/X/indices", f["mod/rna/X/indices"][()].astype(np.int64)), still_read=True)
broken_screen("pointers_64_bit", ["S4.3"], "Column pointers stored as 64-bit integers beside 32-bit cell positions: mixed widths, which scipy rejects in some operations",
              lambda f: replace(f, "mod/rna/X/indptr", f["mod/rna/X/indptr"][()].astype(np.int64)), still_read=True)


def both_wide(f):
    for name in ("indices", "indptr"):
        replace(f, f"mod/rna/X/{name}", f[f"mod/rna/X/{name}"][()].astype(np.int64))


broken_screen("both_index_arrays_64_bit", ["S4.2", "S4.3"], "Both index arrays stored as 64-bit integers in a file small enough for 32-bit. "
              "The same file is what a file past the limit must look like: with the limit lowered for the test, the validators pass it",
              both_wide, still_read=True)
cases[-1]["conforming_with_index_limit"] = 5   # the rna modality has 6 cells
broken_screen("shared_column_named_like_a_modality", ["S5.1"], "A shared cell column named `rna_quality`",
              lambda f: add_column(f, "obs", "rna_quality", np.arange(6, dtype=np.float64)))
broken_screen("covariate_name_with_a_dash", ["S5.3"], "A cell column named `donor-age`, which has a character outside the allowed ones",
              lambda f: add_column(f, "obs", "donor-age", np.arange(6, dtype=np.float64)))


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
broken_screen("grna_group_stored", ["S8.5"], "The gRNA table stores `grna_group`, which is an analysis's choice",
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
def types_changed(change):
    def apply(f):
        types = types_of(f)
        change(types, vectors_of(f)[1])
        add_text_column(f, "grna_type", types)
    return apply


broken_screen("unknown_grna_type", ["S8.15"], "A `grna_type` is not one of the four kinds",
              types_changed(lambda types, targets: types.__setitem__(0, "positive control")), source=high_py)
broken_screen("control_type_with_a_target", ["S8.16"], "A gRNA with a real target is marked `safe-targeting`",
              types_changed(lambda types, targets: types.__setitem__(targets.index(next(t for t in targets if t != "non-targeting")), "safe-targeting")), source=high_py)
broken_screen("targeting_type_without_a_target", ["S8.16"], "A gRNA whose target is `non-targeting` is marked `targeting`",
              types_changed(lambda types, targets: types.__setitem__(targets.index("non-targeting"), "targeting")), source=high_py)
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
broken_screen("positions_without_genome_build", ["S8.14"], "The gRNA table has positions and the file has no `genome_build`",
              lambda f: f["uns"].__delitem__("genome_build"))
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
broken_screen("no_format_entries", ["S10.1"], "The file has neither `perturb_seq_format` nor `dataset_id`",
              lambda f: [f["uns"].__delitem__(name) for name in ("perturb_seq_format", "dataset_id")])


def set_scalar(f, name, value):
    del f[f"uns/{name}"]
    f["uns"].create_dataset(name, data=value, dtype=text if isinstance(value, str) else None)
    encoding = "string" if isinstance(value, str) else "numeric-scalar"
    f[f"uns/{name}"].attrs.update({"encoding-type": encoding, "encoding-version": "0.2.0"})


broken_screen("title_as_a_number", ["S10.5"], "The optional title is a number",
              lambda f: f["uns"].create_dataset("title", data=np.int32(7)) if "title" not in f["uns"] else set_scalar(f, "title", np.int32(7)))
broken_screen("unknown_perturbation_type", ["S10.6"], "The optional perturbation type is not one of the listed kinds",
              lambda f: set_scalar(f, "perturbation_type", "CRISPRx"))
broken_screen("format_from_the_future", ["S10.3"], "The data file's format version is 9.0", lambda f: set_scalar(f, "perturb_seq_format", "9.0"))
broken_screen("format_as_a_number", ["S10.1"], "The data file's format is the integer 0, not a version", lambda f: set_scalar(f, "perturb_seq_format", np.int32(0)))

(ROOT / "cases.json").write_text(json.dumps(cases, indent=2) + "\n")
conforming = sum(1 for c in cases if c["conforms"])
still_read = sum(1 for c in cases if "reads_as" in c)
print(f"{len(cases)} cases: {conforming} conform, {len(cases) - conforming} do not ({still_read} of those a reader still reads); "
      f"rules covered: {len({r for c in cases for r in c['expect']})}")
