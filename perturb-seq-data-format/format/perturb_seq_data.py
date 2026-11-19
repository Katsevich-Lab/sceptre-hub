"""The Perturb-seq data file format (draft, format 0.1): Python validator, and a
temporary writer.

SPEC.md states the rules and RATIONALE.md the reasons. Section numbers below
(S4, S10, ...) refer to both. The structure of a data file is read from
spec/screen.json, which the R code reads too.

    validate_screen()   checks a data file against SPEC.md, building no object
    check_screen()      says whether software must refuse a file it has read
    write_screen()      MuData -> data file; `mudata.read_h5mu` reads it back

Nothing here knows about analyses or about any analysis tool.

Names. A modality is named for what was measured (``rna``), and the gRNA
modality is always ``grna`` (S2).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# The format version, "MAJOR.MINOR" (S10.3). A minor increase adds things an older
# reader can ignore; a major increase is one an older reader must refuse.
FORMAT_VERSION = "0.1"


SPEC_DIR = Path(__file__).resolve().parent.parent / "spec"


SCREEN_SPEC = json.loads((SPEC_DIR / "screen.json").read_text())


# The optional entries of a data file's metadata (S10.5): a name for people, and
# what is known about the screen at import and cannot be worked out later.
METADATA = tuple(SCREEN_SPEC["metadata"]["entries"])


OPTIONAL_METADATA = tuple(k for k, v in SCREEN_SPEC["metadata"]["entries"].items() if not v.get("required"))


VERSION_ENTRY, ID_ENTRY = "perturb_seq_format", "dataset_id"


GRNA = "grna"


# the columns of the gRNA table that hold a target's position (S8.2)
GRNA_COORDINATES = {"chr": "target_chr", "start": "target_start", "end": "target_end"}


# ---------------------------------------------------------------- the data file


def _is_text(column) -> bool:
    return not isinstance(column.dtype, pd.CategoricalDtype) and (
        column.dtype == object or pd.api.types.is_string_dtype(column.dtype)
    )


def _as_file_frame(frame: pd.DataFrame, covariates: bool) -> pd.DataFrame:
    """Column types of the cell and feature tables (S7).
    covariates (cell tables): text and true/false become categoricals
    feature tables:           text stays text, as a categorical only where it has
                              missing values (R's MuData cannot store a missing string)
    both:                     an integer column with missing values becomes doubles"""
    frame = frame.copy()
    frame.index = pd.Index(np.asarray(frame.index, dtype=object))
    for name in frame.columns:
        column = frame[name]
        if pd.api.types.is_bool_dtype(column.dtype):
            labels = column.map({True: "TRUE", False: "FALSE"})
            frame[name] = pd.Categorical(labels, categories=["FALSE", "TRUE"])
        elif _is_text(column):
            if covariates or column.isna().any():
                frame[name] = pd.Categorical(column)
            else:
                frame[name] = column.astype(object)
        elif pd.api.types.is_integer_dtype(column.dtype):
            frame[name] = column.astype("float64") if column.isna().any() else column.astype("int32")
    return frame


def vector_faults(grna_ids, vectors, targets) -> list[tuple[str, str]]:
    """Several gRNAs may be delivered on one vector, named by the optional column
    `vector_id` of the gRNA table (S8.6). A vector has one target, and a vector is
    not named like a gRNA on another vector. Returns the rules broken."""
    faults = []
    seen = {}
    if any(seen.setdefault(v, t) != t for v, t in zip(vectors, targets, strict=True)):
        faults.append(("S8.7", "gRNAs on one vector have different targets; screens that combine targets on a vector are not supported."))
    vector_of = dict(zip(grna_ids, vectors, strict=True))
    if any(v in vector_of and vector_of[v] != v for v in vectors):
        faults.append(("S8.8", "A vector is named like a gRNA on another vector."))
    return faults


GRNA_TYPES = tuple(next(c for c in SCREEN_SPEC["grna_table"] if c["name"] == "grna_type")["values"])


def type_faults(types, targets, placed) -> list[tuple[str, str]]:
    """The optional column `grna_type` of the gRNA table says what kind of gRNA each
    is (S8.15), and agrees with `grna_target`: a control of any kind has the target
    `non-targeting` and no position, and a targeting gRNA has another target
    (S8.16). `placed` says which gRNAs have a target position. Returns the rules broken."""
    if any(t not in GRNA_TYPES for t in types):
        return [("S8.15", f"`grna_type` must be one of {', '.join(GRNA_TYPES)}.")]
    control = [t != "targeting" for t in types]
    if any(c != (target == "non-targeting") or (c and has) for c, target, has in zip(control, targets, placed, strict=True)):
        return [("S8.16", "`grna_type` and `grna_target` disagree: a control of any kind has the target `non-targeting` and no position, and a targeting gRNA has another target.")]
    return []


def index_type(n_stored, n_cells, n_features, index_limit=None) -> str:
    """The one type of both index arrays of a count matrix (S4.2, S4.3): 32-bit
    while the stored values, the cells and the features each number at most
    2,147,483,647, and 64-bit beyond. It is the type scipy chooses for a matrix of
    that size, so a matrix scipy built already has it. `index_limit` replaces the
    limit, for tests only: no test file is large enough to reach the real one."""
    rule = SCREEN_SPEC["counts"]["index"]
    limit = rule["narrow_up_to"] if index_limit is None else index_limit
    return rule["narrow"] if max(int(n_stored), int(n_cells), int(n_features)) <= limit else rule["wide"]


def write_screen(mdata, file, index_limit=None) -> None:
    """Write a MuData as a data file. `uns` holds the entries of S10 that are
    known, with `dataset_id` if the dataset already has one. Every modality is written
    under its own name; `grna` is the gRNA modality."""
    import anndata as ad
    import mudata as md
    from scipy import sparse

    file = Path(file)
    if file.exists():
        raise ValueError(f"{file} exists; a data file is written once.")
    stamp = {key: mdata.uns[key] for key in METADATA if key in mdata.uns}
    grna_table = mdata.mod[GRNA].var
    if not stamp.get("genome_build") and any(grna_table[c].notna().any() for c in GRNA_COORDINATES.values() if c in grna_table):
        raise ValueError("The gRNA table has positions: uns['genome_build'] must say which genome build they are in.")
    names = list(mdata.mod)
    if GRNA not in names or len(names) < 2:
        raise ValueError(f"A data file needs the modality '{GRNA}' and at least one other.")
    if not all(re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", n) for n in names):
        raise ValueError("Modality names must start with a letter and use letters, digits and underscores.")
    cells = pd.Index(np.asarray(mdata[GRNA].obs_names, dtype=object))
    features = [f for n in names for f in mdata[n].var_names]
    if not cells.is_unique:
        raise ValueError("Cell identifiers must be unique.")
    if len(set(features)) != len(features):
        raise ValueError("Feature identifiers must be unique across the modalities of a file.")
    targets = mdata[GRNA].var.get("grna_target")
    if targets is None or targets.isna().any():
        raise ValueError("Every gRNA needs a `grna_target`.")
    _check_vectors(mdata[GRNA].var, None)
    if "grna_type" in mdata[GRNA].var:
        placed = np.zeros(mdata[GRNA].n_vars, dtype=bool)
        for column in GRNA_COORDINATES.values():
            if column in mdata[GRNA].var:
                placed |= mdata[GRNA].var[column].notna().to_numpy()
        kinds = [None if pd.isna(t) else str(t) for t in mdata[GRNA].var["grna_type"]]
        for _, message in type_faults(kinds, list(map(str, targets)), placed):
            raise ValueError(message)
    modalities = {}
    for name in names:
        source = mdata[name]
        if list(source.obs_names) != list(cells):
            raise ValueError("The modalities must hold the same cells in the same order.")
        # column-compressed cells x features: one feature is one contiguous read (S3)
        counts = sparse.csc_matrix(source.X)
        counts.sum_duplicates()
        counts.eliminate_zeros()
        counts.sort_indices()
        if not np.array_equal(counts.data, np.round(counts.data)):
            raise ValueError("Counts must be integers.")
        # every modality's cell table has its total count and its number of features
        # with a nonzero count; a writer computes them if they are absent (S6)
        obs = source.obs.loc[cells].copy()
        if "n_umis" not in obs:
            obs["n_umis"] = np.asarray(counts.sum(axis=1)).reshape(-1).astype(np.float64)
        if "n_nonzero" not in obs:
            obs["n_nonzero"] = np.bincount(counts.indices, minlength=len(cells)).astype(np.int32)
        # counts as 32-bit integers, and both index arrays in the one type the size
        # requires (S4), set before writing: an array rewritten afterwards would
        # leave its first copy as dead space in the file
        counts.data = counts.data.astype(np.int32)
        width = np.dtype(index_type(counts.nnz, *counts.shape, index_limit))
        counts.indices = counts.indices.astype(width)
        counts.indptr = counts.indptr.astype(width)
        modalities[name] = ad.AnnData(X=counts, obs=_as_file_frame(obs, True), var=_as_file_frame(source.var, False))
    out = md.MuData(modalities)
    shared = _as_file_frame(mdata.obs.loc[cells, [c for c in mdata.obs.columns if ":" not in str(c)]], True)
    for name in shared.columns:
        out.obs[name] = shared[name].array
    # scalar entries at the top level of `uns` (S10)
    out.uns[VERSION_ENTRY] = FORMAT_VERSION
    out.uns[ID_ENTRY] = str(stamp.get(ID_ENTRY) or uuid.uuid4().hex)
    # optional: a title for people, and what is known about the screen (S10.5)
    for key in OPTIONAL_METADATA:
        if stamp.get(key):
            out.uns[key] = str(stamp[key])
    file.parent.mkdir(parents=True, exist_ok=True)
    # plain HDF5 strings: R cannot open a file with anndata's nullable string arrays (S7)
    # noted before writing: mudata converts these columns in the object as it writes
    text_columns = {
        name: {c: np.asarray(out[name].var[c], dtype=object) for c in out[name].var.columns if _is_text(out[name].var[c]) and not out[name].var[c].isna().any()}
        for name in out.mod
    }
    with ad.settings.override(allow_write_nullable_strings=False):
        out.write(file)
    _text_as_text(file, text_columns)


def _text_as_text(file, text_columns: dict) -> None:
    """A text column of a feature table without missing values is stored as text
    (S7.2), as the R writer stores it. mudata turns a text column with repeated
    values into a categorical when it writes; this puts such columns back."""
    import anndata as ad
    import h5py

    with h5py.File(file, "r+") as f:
        for name, columns in text_columns.items():
            group = f[f"mod/{name}/var"]
            for column, values in columns.items():
                if group[column].attrs.get("encoding-type") == "string-array":
                    continue
                del group[column]
                with ad.settings.override(allow_write_nullable_strings=False):
                    ad.io.write_elem(group, column, values)


def format_parts(version):
    """(major, minor) of a format version, or None if it is not one."""
    if isinstance(version, bytes):
        version = version.decode()
    if not isinstance(version, str) or not re.fullmatch(r"[0-9]+\.[0-9]+", version):
        return None
    major, minor = version.split(".")
    return int(major), int(minor)


def known_major(version) -> bool:
    parts = format_parts(version)
    return parts is not None and parts[0] == format_parts(FORMAT_VERSION)[0]


def _check_vectors(grna_table: pd.DataFrame, file) -> None:
    if "vector_id" not in grna_table:
        return
    if grna_table["vector_id"].isna().any():
        raise ValueError("Every gRNA needs a `vector_id` when the column is present.")
    for _, message in vector_faults(list(grna_table.index), list(map(str, grna_table["vector_id"])), list(map(str, grna_table["grna_target"]))):
        raise ValueError(f"{file}: {message}" if file is not None else message)


def check_format(version, what) -> None:
    if format_parts(version) is None:
        raise ValueError(f"{what} has no format version.")
    if not known_major(version):
        raise ValueError(f"{what} is in format {version}; this software knows format {FORMAT_VERSION} and its minor versions.")


def check_screen(mdata) -> None:
    """What software must refuse before using a data file it has read with
    `mudata.read_h5mu`: a major format version it does not know (S10.3), and a
    gRNA table whose vectors break S8.7 or S8.8 (S8.9). Cheap: it looks at the
    metadata and the gRNA table only."""
    version = mdata.uns.get(VERSION_ENTRY)
    check_format(None if version is None else str(np.asarray(version).reshape(-1)[0]), "The data file")
    _check_vectors(mdata[GRNA].var, "The data file")


# ---------------------------------------------------------------- validation

# validate_screen() checks a data file against the numbered rules of SPEC.md. It
# reads the file directly, builds no object, and reports every problem it finds: a
# DataFrame with one row per problem (`rule`, `where`, `message`), with no rows
# when the file conforms. The R library has the same function and reports the same
# rules.

MODALITY_NAME = SCREEN_SPEC["modalities"]["name_pattern"]


class _Problems:
    def __init__(self):
        self.found = set()

    def add(self, rule, where, message):
        self.found.add((rule, str(where), message))

    def all(self) -> pd.DataFrame:
        return pd.DataFrame(sorted(self.found), columns=["rule", "where", "message"])


def _storage_kind(node) -> str:
    """The kind of a column of a cell or feature table (spec/screen.json, storage_kinds)."""
    import h5py

    encoding = _encoding(node)
    if encoding == "categorical":
        return "categorical"
    if encoding == "string-array":
        return "string"
    if encoding == "array" and isinstance(node, h5py.Dataset) and node.dtype == np.float64:
        return "double"
    if encoding == "array" and isinstance(node, h5py.Dataset) and node.dtype == np.int32:
        return "int32"
    return "other"


def _text(values) -> list:
    return [v.decode() if isinstance(v, bytes) else str(v) for v in np.asarray(values).reshape(-1)]


def _encoding(node) -> str | None:
    value = node.attrs.get("encoding-type")
    return value.decode() if isinstance(value, bytes) else value


def _filters(dataset) -> set:
    import h5py

    plist = dataset.id.get_create_plist()
    names = {h5py.h5z.FILTER_DEFLATE: "deflate", h5py.h5z.FILTER_SHUFFLE: "shuffle"}
    return {names.get(plist.get_filter(i)[0], "other") for i in range(plist.get_nfilters())}


def _index_of(table):
    name = table.attrs.get("_index")
    name = name.decode() if isinstance(name, bytes) else name
    return _text(table[name][()]) if name is not None and name in table else None


def validate_screen(file, index_limit=None) -> pd.DataFrame:
    """`index_limit` is for tests only: see `index_type`."""
    import h5py

    log = _Problems()
    try:
        f = h5py.File(file, "r")
    except OSError:
        log.add("S1.1", file, "is not an HDF5 file.")
        return log.all()
    with f:
        if _encoding(f) != SCREEN_SPEC["root"]["encoding-type"]:
            log.add("S1.1", file, "The file is not marked as MuData.")
        if "mod" not in f:
            log.add("S2.1", file, "The file has no modalities.")
            return log.all()
        in_file = sorted(f["mod"])
        if GRNA not in in_file or len(in_file) < SCREEN_SPEC["modalities"]["minimum"]:
            log.add("S2.1", "mod", "A data file needs the modality 'grna' and at least one other.")
        for name in in_file:
            if not re.fullmatch(MODALITY_NAME.strip("^$"), name):
                log.add("S2.2", f"mod/{name}", "is not a valid modality name.")

        # file-level metadata
        where = SCREEN_SPEC["metadata"]["path"]
        entries = SCREEN_SPEC["metadata"]["entries"]
        record = f.get(where)

        def scalar_of(name, kind):
            kinds = "iu" if kind == "integer" else "OSU"
            return record is not None and name in record and isinstance(record[name], h5py.Dataset) and record[name].shape == () and record[name].dtype.kind in kinds

        for name, entry in entries.items():
            if not entry["required"] and record is not None and name in record and not scalar_of(name, entry["kind"]):
                log.add("S10.5", f"{where}/{name}", "must be a scalar text entry.")
        if not all(scalar_of(name, entry["kind"]) for name, entry in entries.items() if entry["required"]):
            log.add("S10.1", where, f"needs the scalar text entries `{VERSION_ENTRY}` and `{ID_ENTRY}`.")
        elif format_parts(_text(record[VERSION_ENTRY][()])[0]) is None:
            log.add("S10.1", f"{where}/{VERSION_ENTRY}", "is not a version of the form MAJOR.MINOR.")
        else:
            if not known_major(_text(record[VERSION_ENTRY][()])[0]):
                log.add("S10.3", f"{where}/{VERSION_ENTRY}", "is a format this validator does not know.")
            if not _text(record[ID_ENTRY][()])[0]:
                log.add("S10.1", f"{where}/{ID_ENTRY}", "is empty.")
        # optional entries with a fixed vocabulary (S10.6)
        for name in OPTIONAL_METADATA:
            allowed = entries[name].get("values")
            if allowed and scalar_of(name, "string") and _text(record[name][()])[0] not in allowed:
                log.add("S10.6", f"{where}/{name}", f"must be one of {', '.join(allowed)}.")

        # cells: the same, in the same order, in every modality and in the shared table
        cells = _index_of(f[f"mod/{GRNA if GRNA in in_file else in_file[0]}/obs"])
        if cells is None:
            log.add("S2.3", "mod", "A modality has no cell identifiers.")
            return log.all()
        if len(set(cells)) != len(cells):
            log.add("S9.1", "obs", "Cell identifiers are not unique.")
        for name in in_file:
            if _index_of(f[f"mod/{name}/obs"]) != cells:
                log.add("S2.3", f"mod/{name}/obs", "does not hold the same cells in the same order as 'grna'.")
        if "obs" in f and _index_of(f["obs"]) != cells:
            log.add("S2.3", "obs", "The shared cell table does not hold the same cells in the same order.")
        features = {name: _index_of(f[f"mod/{name}/var"]) or [] for name in in_file}
        every = [x for name in in_file for x in features[name]]
        if len(set(every)) != len(every):
            log.add("S9.2", "mod", "Feature identifiers are not unique across the modalities.")
        if any("\n" in x for x in cells + every):
            log.add("S9.4", "mod", "An identifier contains a line feed.")

        # the count matrix of every modality (S3, S4, S11)
        expect, rules = SCREEN_SPEC["counts"]["arrays"], SCREEN_SPEC["counts"]["rules"]
        for name in in_file:
            x = f"mod/{name}/X"
            group = f.get(x)
            shape = None if group is None else group.attrs.get("shape")
            if group is None or _encoding(group) != SCREEN_SPEC["counts"]["encoding-type"] or shape is None or [int(v) for v in shape] != [len(cells), len(features[name])]:
                log.add("S3.1", x, "is not a column-compressed (csc_matrix) matrix of cells x features.")
                continue
            # both index arrays in the one type this file's size requires (S4.2, S4.3)
            index = index_type(group["data"].shape[0], len(cells), len(features[name]), index_limit)
            for array, kind in expect.items():
                kind, why = (index, " (the index type a file of this size requires)") if kind == "index" else (kind, "")
                if group[array].dtype != np.dtype(kind):
                    # not conforming; a reader still accepts any whole-number type here (S4.4)
                    log.add(rules[array], f"{x}/{array}", f"is not stored as {kind.removeprefix('int')}-bit signed integers{why}.")
                if _filters(group[array]) - set(SCREEN_SPEC["counts"]["filters"]):
                    log.add("S11.1", f"{x}/{array}", "uses a filter other than gzip.")
            data, indices, indptr = (group[k][()] for k in ("data", "indices", "indptr"))
            if len(indptr) != len(features[name]) + 1 or indptr[0] != 0 or np.any(np.diff(indptr.astype(np.int64)) < 0) or indptr[-1] != len(data) or len(indices) != len(data):
                log.add("S3.2", x, "The column pointers do not describe the stored values.")
                continue
            if len(indices):
                step = np.diff(indices.astype(np.int64))
                starts = indptr[(indptr > 0) & (indptr < len(indices))].astype(np.int64)   # where a new feature begins
                within = np.delete(step, starts - 1)
                if np.any((indices < 0) | (indices >= len(cells))) or np.any(within <= 0):
                    log.add("S3.3", x, "Within a feature, cell positions must be in range and strictly increasing.")
                if np.any(data <= 0) or np.any(data != np.round(data)):
                    log.add("S3.4", x, "Stored counts must be positive integers.")

        # column kinds (S7); returns each column's kind
        def check_table(table, cell_table):
            allowed = SCREEN_SPEC["cell_tables" if cell_table else "feature_tables"]["column_kinds"]
            group = f[table]
            index = group.attrs.get("_index")
            index = index.decode() if isinstance(index, bytes) else index
            kinds = {}
            for name in group:
                if name == index:
                    continue
                path, node = f"{table}/{name}", group[name]
                kind = kinds[name] = _storage_kind(node)
                if kind not in allowed:
                    log.add("S7.1" if cell_table else "S7.2", path, f"is not one of the column kinds allowed here ({', '.join(allowed)}).")
                elif kind == "categorical":
                    categories, codes = node["categories"][()], node["codes"][()]
                    texts = _text(categories)
                    if categories.dtype.kind not in "OSU" or len(set(texts)) != len(texts) or np.any((codes < -1) | (codes >= len(texts))):
                        log.add("S7.3", path, "is not a well-formed categorical.")
                elif kind == "int32":
                    if np.any(node[()] == -2147483648):
                        log.add("S7.4", path, "An integer column holds a missing value; it must be stored as doubles.")
                if cell_table and not re.fullmatch(SCREEN_SPEC["cell_tables"]["column_name_pattern"].strip("^$"), name):
                    log.add("S5.3", path, "is not a valid column name (letters, digits, underscores and dots, starting with a letter).")
            return kinds

        if "obs" in f:
            shared = check_table("obs", True)
            for name in in_file:
                for column in shared:
                    if column.startswith(f"{name}_"):
                        log.add("S5.1", f"obs/{column}", "A shared column must not start with a modality's name and an underscore.")
        for name in in_file:
            table = f"mod/{name}/obs"
            kinds = check_table(table, True)
            # the columns every modality's cell table has (S6)
            for column in SCREEN_SPEC["cell_tables"]["every_modality"]:
                if column["name"] not in kinds:
                    if column["required"]:
                        log.add("S6.1", table, f"has no column `{column['name']}`.")
                elif kinds[column["name"]] not in column["kinds"]:
                    log.add("S6.2", f"{table}/{column['name']}", "must be a number.")
            check_table(f"mod/{name}/var", False)

        # the gRNA table (S8)
        if GRNA in in_file:
            table = f"mod/{GRNA}/var"
            group = f[table]

            def text(name):
                node = group[name]
                if _encoding(node) == "string-array":
                    return _text(node[()])
                categories = _text(node["categories"][()])
                return [None if c < 0 else categories[c] for c in node["codes"][()]]

            for column in SCREEN_SPEC["grna_table_forbidden"]:
                if column["name"] in group:
                    log.add(column["rule"], f"{table}/{column['name']}", "must not be stored.")
            for column in SCREEN_SPEC["grna_table"]:
                path = f"{table}/{column['name']}"
                if column["name"] not in group:
                    if column["required"]:
                        log.add(column["rule"], table, f"The gRNA table has no `{column['name']}`.")
                    continue
                if _storage_kind(group[column["name"]]) not in column["kinds"]:
                    what = "text" if "string" in column["kinds"] else "a number"
                    log.add(column["rule"], path, f"must be {what}" + (" with no missing values." if column["required"] else "."))
                elif column["required"] or column.get("no_missing"):
                    if any(t is None or t == "" for t in text(column["name"])):
                        log.add(column["rule"], path, "must be text with no missing values.")
            # target positions: 1-based, both ends included; all three or none (S8.10, S8.11)
            n_grnas = len(features[GRNA])

            def position(name, as_text):
                if name not in group:
                    return [None] * n_grnas if as_text else np.full(n_grnas, np.nan)
                kind = _storage_kind(group[name])
                if as_text:
                    return text(name) if kind in ("string", "categorical") else None
                return group[name][()].astype(float) if kind in ("double", "int32") else None

            chromosome, start, end = position("target_chr", True), position("target_start", False), position("target_end", False)
            if chromosome is not None and start is not None and end is not None:
                given = np.array([c is not None for c in chromosome]).astype(int) + (~np.isnan(start)).astype(int) + (~np.isnan(end)).astype(int)
                if np.any((given == 1) | (given == 2)):
                    log.add("S8.11", table, "A gRNA has some of `target_chr`, `target_start`, `target_end` and not all three.")
                if np.any(given > 0) and where in f and "genome_build" not in f[where]:
                    log.add("S8.14", table, "The gRNA table has positions and the file has no `genome_build`.")
                both = ~np.isnan(start) & ~np.isnan(end)
                a, b = start[both], end[both]
                if np.any((a != np.round(a)) | (b != np.round(b)) | (a < 1) | (b < a)):
                    log.add("S8.10", table, "Target positions must be whole numbers, 1-based with both ends included: `target_start` at least 1 and `target_end` not before it.")
            # gRNAs delivered on one vector (S8.6 to S8.8)
            if "vector_id" in group and "grna_target" in group:
                kinds = {"string", "categorical"}
                if _storage_kind(group["vector_id"]) in kinds and _storage_kind(group["grna_target"]) in kinds:
                    vectors, targets = text("vector_id"), text("grna_target")
                    if None not in vectors and None not in targets:
                        for rule, message in vector_faults(features[GRNA], vectors, targets):
                            log.add(rule, f"{table}/vector_id", message)
            # the kind of each gRNA, where recorded (S8.15, S8.16)
            if "grna_type" in group and "grna_target" in group:
                kinds = {"string", "categorical"}
                if _storage_kind(group["grna_type"]) in kinds and _storage_kind(group["grna_target"]) in kinds and None not in text("grna_target"):
                    placed = np.zeros(n_grnas, dtype=bool)
                    if chromosome is not None:
                        placed |= np.array([c is not None for c in chromosome])
                    for values in (start, end):
                        if values is not None:
                            placed |= ~np.isnan(values)
                    for rule, message in type_faults(text("grna_type"), text("grna_target"), placed):
                        log.add(rule, f"{table}/grna_type", message)
    return log.all()
