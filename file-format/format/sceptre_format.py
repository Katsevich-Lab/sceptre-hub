"""The sceptre file format (draft, format 0.1): Python reader and writer.

SPEC.md states the rules and RATIONALE.md the reasons. Section numbers below
(S4, T3, ...) refer to both.

Two layers, as in sceptre_format.R:
  1. Files <-> standard Python containers, with no use of pysceptre:
       write_screen(), read_screen()   data file <-> MuData
       init_analysis(), write_stage(), read_analysis(), clear_stage()
  2. Checking files against SPEC.md, without building any object:
       validate_screen(), validate_analysis()

Nothing in this file uses pysceptre.

The tables of an analysis, the stages and the shape of the JSON records are not
written here: they are read from spec/ (tables.json, screen.json and two JSON
Schemas), which the R library reads too.

Names. The data file names a modality for what was measured (``rna``), and the
gRNA modality is always ``grna`` (S2). An analysis says which modality is its
response (A1).
"""

from __future__ import annotations

import hashlib
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
TABLE_SPEC = json.loads((SPEC_DIR / "tables.json").read_text())
SCREEN_SPEC = json.loads((SPEC_DIR / "screen.json").read_text())
# The optional entries of a data file's metadata (S10.5): a name for people, and
# what is known about the screen at import and cannot be worked out later.
OPTIONAL_METADATA = tuple(k for k, v in SCREEN_SPEC["metadata"]["entries"].items() if not v.get("required"))
GRNA = "grna"
DEFAULT_RESPONSE_MODALITY = "rna"
# the gRNA table's target coordinates: name in the package -> name in the file (S8)
GRNA_COORDINATES = {"chr": "target_chr", "start": "target_start", "end": "target_end"}
STAGES = tuple(TABLE_SPEC["stages"])
# calibration, power and discovery do not depend on one another (A4)
STAGE_RANK = dict(TABLE_SPEC["stage_rank"])
STAGE_FUNCTIONS = {
    "parameters": "set_analysis_parameters",
    "assignment": "assign_grnas",
    "qc": "run_qc",
    "calibration": "run_calibration_check",
    "power": "run_power_check",
    "discovery": "run_discovery_analysis",
}
# The shared list of response-fitting methods (T6). A method fitted in one
# language only is never reused by the other.
FIT_METHODS = {
    "batched": ("R", "Python"),
    "sceptre": ("R", "Python"),
    "pysceptre": ("R", "Python"),
    "glmGamPoi": ("R",),
}
ROW_GROUP_ROWS = int(TABLE_SPEC["max_rows_per_row_group"])

# ---------------------------------------------------------------- small helpers


def sha256_ids(ids) -> str:
    return hashlib.sha256("\n".join(map(str, ids)).encode()).hexdigest()


def _atomically(path, write) -> None:
    """Each file is written under a temporary name and renamed into place (A6)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    write(temporary)
    os.replace(temporary, path)


def _write_json(value, path) -> None:
    _atomically(path, lambda t: Path(t).write_text(json.dumps(value, indent=2) + "\n"))


def _read_json(path):
    return json.loads(Path(path).read_text())


def _byte_order(table: pd.DataFrame, key) -> list[int]:
    """Row order by the key columns, each compared as UTF-8 bytes, missing last.
    The sort is stable, and done by numpy on byte strings, not row by row."""
    columns = []
    for name in reversed(list(key)):   # lexsort takes the most significant key last
        values = table[name]
        missing = values.isna().to_numpy()
        encoded = np.array([b"" if m else str(v).encode() for v, m in zip(values, missing)], dtype="S")
        columns += [encoded, missing]
    if not columns or len(table) == 0:
        return list(range(len(table)))
    return np.lexsort(columns).tolist()


# ---------------------------------------------------------------- tables (T0)


def write_table(table: pd.DataFrame, path, key, metadata=None) -> int:
    """Rows sorted by the key in byte order; row groups of at most a million rows.
    `metadata` is text by key, stored in the file's own key-value metadata."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    table = table.reset_index(drop=True).copy()
    for name in table.columns:
        column = table[name]
        if isinstance(column.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(column.dtype):
            table[name] = column.astype(object).where(column.notna(), None)
        elif pd.api.types.is_bool_dtype(column.dtype):
            table[name] = column.astype("boolean")
        elif pd.api.types.is_integer_dtype(column.dtype):
            table[name] = column.astype("int32")
    table = table.iloc[_byte_order(table, key)].reset_index(drop=True)
    arrays, fields = [], []
    for name in table.columns:
        column = table[name]
        if column.dtype == object and len(column) and isinstance(column.iloc[0], (list, np.ndarray)):
            array = pa.array([None if v is None else list(map(float, v)) for v in column], type=pa.list_(pa.float64()))
        elif column.dtype == object:
            array = pa.array(column.tolist(), type=pa.string())
        else:
            array = pa.Array.from_pandas(column)
        arrays.append(array)
        fields.append(pa.field(str(name), array.type))
    out = pa.Table.from_arrays(arrays, schema=pa.schema(fields, metadata=metadata))
    _atomically(path, lambda t: pq.write_table(out, t, row_group_size=ROW_GROUP_ROWS))
    return len(table)


def _table_files(path) -> list[Path]:
    if isinstance(path, (list, tuple)):
        return [Path(f) for f in path]
    path = Path(path)
    return sorted(path.glob("*.parquet")) if path.is_dir() else [path]


def table_rows(path) -> int:
    import pyarrow.parquet as pq

    return sum(pq.ParquetFile(f).metadata.num_rows for f in _table_files(path))


def read_table(path, columns=None) -> pd.DataFrame:
    """One .parquet file or a folder of part files; only the known columns where
    `columns` is given; identifiers as text whether stored plain or as a dictionary."""
    import pyarrow.parquet as pq

    parts = []
    for file in _table_files(path):
        part = pq.read_table(file)
        if columns is not None:
            part = part.select([c for c in columns if c in part.column_names])
        parts.append(part.to_pandas())
    table = pd.concat(parts, ignore_index=True) if len(parts) > 1 else parts[0]
    for name in table.columns:
        column = table[name]
        if isinstance(column.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(column.dtype):
            table[name] = column.astype(object).where(column.notna(), None)
    return table


def read_table_from_arrow(part) -> pd.DataFrame:
    """One Parquet part as `read_table` returns it: identifiers as plain text."""
    table = part.to_pandas()
    for name in table.columns:
        column = table[name]
        if isinstance(column.dtype, pd.CategoricalDtype) or pd.api.types.is_string_dtype(column.dtype):
            table[name] = column.astype(object).where(column.notna(), None)
    return table


def unit_key(table) -> list[str]:
    return [c for c in ("grna_target", "grna_id", "response_id") if c in table.columns]


# Every table of an analysis, from spec/tables.json: the stage that writes it,
# its key, and the columns a reader takes (None: all of them).
TABLE_DEFINITIONS = {t["name"]: t for t in TABLE_SPEC["tables"]}
TABLES = {
    name: (t["stage"], list(t["key"]), [c["name"] for c in t["columns"]] if t["extra_columns"] == "ignored" else None)
    for name, t in TABLE_DEFINITIONS.items()
}


def _table_key(name, table):
    """A key column that is not required is part of the key when it is present."""
    return [k for k in TABLES[name][1] if k in table.columns]


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


def _from_file_frame(frame: pd.DataFrame, covariates: bool) -> pd.DataFrame:
    frame = frame.copy()
    for name in frame.columns:
        column = frame[name]
        if not covariates and (isinstance(column.dtype, pd.CategoricalDtype) or _is_text(column)):
            frame[name] = column.astype(object).where(column.notna(), None)
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


def response_modality_of(modalities, requested=None, what="The data file") -> str:
    """The modality an analysis tests, among those of a data file: the one asked
    for, or the only one besides ``grna`` (S2)."""
    modalities = list(modalities)
    if GRNA not in modalities:
        raise ValueError(f"{what} has no modality named '{GRNA}'.")
    candidates = [m for m in modalities if m != GRNA]
    if requested is None:
        if len(candidates) == 1:
            return candidates[0]
        raise ValueError(f"{what} has the modalities {', '.join(candidates)} besides '{GRNA}'; say which one is the response.")
    if requested not in candidates:
        raise ValueError(f"{what} has no modality named '{requested}'.")
    return requested


def response_of(mdata) -> str:
    """The name of the response modality of a MuData that `read_screen` returned."""
    return mdata.uns["sceptre"]["response_modality"]


def covariate_table(mdata) -> pd.DataFrame:
    """The model's covariate table (S5). Facts about the cell are in the shared cell
    table, `mdata.obs`; measurements of a modality are in that modality's `obs`. The
    table is all of them, a modality's columns prefixed with its name and an
    underscore (`rna_n_umis`, `grna_n_umis`). Only the response modality and `grna`
    contribute."""
    used = (response_of(mdata), GRNA)
    cells = mdata[used[0]].obs_names
    frames = [mdata[m].obs.loc[cells].add_prefix(f"{m}_") for m in used]
    shared = mdata.obs.loc[cells, [c for c in mdata.obs.columns if ":" not in str(c)]]
    frames.append(shared)
    return pd.concat([f for f in frames if f.shape[1]] or [shared], axis=1)


def split_covariates(covariates: pd.DataFrame, modalities) -> dict:
    """The inverse of `covariate_table`."""
    out, taken = {}, []
    for m in modalities:
        names = [c for c in covariates.columns if str(c).startswith(f"{m}_")]
        out[m] = covariates[names].rename(columns=lambda c, m=m: c[len(m) + 1 :])
        taken += names
    out["shared"] = covariates.drop(columns=taken)
    return out


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
    """Write a MuData (as `read_screen` returns it) as a data file. `uns["sceptre"]`
    holds `moi`, and `id` if the dataset already has one. Every modality is written
    under its own name; `grna` is the gRNA modality."""
    import anndata as ad
    import mudata as md
    from scipy import sparse

    file = Path(file)
    if file.exists():
        raise ValueError(f"{file} exists; a data file is written once.")
    stamp = dict(mdata.uns["sceptre"])
    if stamp.get("moi") not in ("high", "low"):
        raise ValueError("uns['sceptre']['moi'] must be 'high' or 'low'.")
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
    # a native record with scalar entries (S10)
    out.uns["sceptre"] = {
        "format": FORMAT_VERSION,
        "moi": stamp["moi"],
        "id": stamp.get("id") or uuid.uuid4().hex,
    }
    # optional: a label for people, and where the data came from (S10.5)
    for key in OPTIONAL_METADATA:
        if stamp.get(key):
            out.uns["sceptre"][key] = str(stamp[key])
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
        raise ValueError(f"{what} has no sceptre format version.")
    if not known_major(version):
        raise ValueError(f"{what} is in format {version}; this reader knows format {FORMAT_VERSION} and its minor versions.")


def _scalar(value):
    value = np.asarray(value).reshape(-1)[0]
    return value.item() if isinstance(value, np.generic) else value


def _same_index_width(counts) -> None:
    """A reader accepts any whole-number type for the index arrays (S4.4). A conforming file already holds both
    index arrays in one type, and this does nothing to it. A file from another
    writer may hold them in two, which scipy rejects in some operations; both are
    then brought to the type the size requires."""
    width = np.dtype(index_type(counts.nnz, *counts.shape))
    if counts.indptr.dtype != width or counts.indices.dtype != width:
        counts.indptr = counts.indptr.astype(width, copy=False)
        counts.indices = counts.indices.astype(width, copy=False)


def read_screen(file, response_modality=None):
    """A data file as a MuData: its modalities under their own names (cells x
    features), cell facts in ``mdata.obs``, modality measurements in each modality's
    ``obs``, the gRNA table in ``mdata["grna"].var``, and ``mdata.uns["sceptre"]`` =
    format, moi, id, response_modality.

    `response_modality` says which modality an analysis will test; it may be left
    out when the file has one modality besides ``grna``. Modalities other than that
    one and ``grna`` are carried along and otherwise ignored."""
    import mudata as md

    mdata = md.read_h5mu(file)
    stamp = {k: _scalar(v) for k, v in dict(mdata.uns.get("sceptre", {})).items()}
    check_format(stamp.get("format"), file)
    main = response_modality_of(mdata.mod, response_modality, file)
    if list(mdata[main].obs_names) != list(mdata[GRNA].obs_names):
        raise ValueError("The modalities must hold the same cells in the same order.")
    # every modality, used by an analysis or not, so that the object can be written again
    for name in mdata.mod:
        mdata[name].obs = _from_file_frame(mdata[name].obs, True)
        mdata[name].var = _from_file_frame(mdata[name].var, False)
        _same_index_width(mdata[name].X)
    mdata.obs = _from_file_frame(mdata.obs, True)
    # a reader refuses a file whose vectors break S8.7 or S8.8 (S8.9)
    _check_vectors(mdata[GRNA].var, file)
    mdata.uns["sceptre"] = {
        "format": str(stamp["format"]),
        "moi": str(stamp["moi"]),
        "id": str(stamp["id"]),
        "response_modality": main,
        **{key: str(stamp[key]) for key in OPTIONAL_METADATA if key in stamp},
    }
    return mdata


def screen_link(mdata, file=None, directory=None) -> dict:
    """What an analysis records about its data file (A1)."""
    response = mdata[response_of(mdata)]
    link = {
        "id": mdata.uns["sceptre"]["id"],
        "n_cells": int(response.n_obs),
        "n_responses": int(response.n_vars),
        "n_grnas": int(mdata[GRNA].n_vars),
        "cell_ids_sha256": sha256_ids(response.obs_names),
        "response_ids_sha256": sha256_ids(response.var_names),
        "grna_ids_sha256": sha256_ids(mdata[GRNA].var_names),
    }
    # for people: which dataset this is; not part of the identity a reader checks
    if mdata.uns["sceptre"].get("name"):
        link["name"] = mdata.uns["sceptre"]["name"]
    # a hint only, and only when the data file is near the analysis (A1.4)
    if file is not None and nearby_path(file, directory) is not None:
        link["path"] = nearby_path(file, directory)
    return link


def nearby_path(file, directory):
    """The data file's path relative to an analysis folder, if the file is in the
    folder's parent or grandparent (the `analyses/<name>/` convention) or below one
    of them; otherwise None. A far-away file would give a long chain of `..` that
    names directories and helps nobody."""
    path = Path(os.path.relpath(Path(file).resolve(), Path(directory).resolve())).as_posix()
    return None if path.split("/").count("..") > 2 else path


# ---------------------------------------------------------------- the analysis folder


def init_analysis(directory, link: dict, response_modality: str) -> Path:
    """Start an analysis folder: which data file it analyses and which of its
    modalities is the response (A1)."""
    directory = Path(directory)
    if directory.exists() and any(directory.iterdir()):
        raise ValueError(f"{directory} is not empty.")
    _write_json(
        {"format": FORMAT_VERSION, "response_modality": response_modality, "screen": link},
        directory / "analysis.json",
    )
    return directory


def analysis_response_modality(directory) -> str:
    """The modality an analysis folder names as its response."""
    name = _read_json(Path(directory) / "analysis.json").get("response_modality")
    if name is None:
        raise ValueError(f"{directory} does not say which modality is its response.")
    return name


def _stage_file(directory, stage) -> Path:
    """An analysis folder has one folder per stage: the stage's record and its
    tables (A2). A table's name is its path without `.parquet`: "qc/cells"."""
    return Path(directory) / stage / "record.json"


def table_path(directory, name) -> Path:
    return Path(directory) / f"{name}.parquet"


def stages_run(directory) -> list[str]:
    return [s for s in STAGES if _stage_file(directory, s).exists()]


def write_stage(directory, stage, record: dict, tables=None, provenance=None) -> None:
    """Write one stage into its folder: its tables first, then its record, which
    lists them (A2, A5)."""
    directory = Path(directory)
    if not (directory / "analysis.json").exists():
        raise ValueError(f"{directory} is not an analysis folder.")
    missing = [s for s in STAGES if STAGE_RANK[s] < STAGE_RANK[stage] and s not in stages_run(directory)]
    if missing:
        raise ValueError(f"Stage '{stage}' needs the stage(s) {', '.join(missing)} first.")
    if _stage_file(directory, stage).exists():
        raise ValueError(f"Stage '{stage}' is already written; clear it first.")
    # a folder without a record is a stage that did not finish: nothing of it is kept (A4.1)
    shutil.rmtree(directory / stage, ignore_errors=True)
    files = []
    for name, table in (tables or {}).items():
        if name not in TABLES or TABLES[name][0] != stage:
            raise ValueError(f"Stage '{stage}' does not write the table {name}.")
        # listed by its name within the stage's folder
        files.append({"path": f"{name.split('/')[-1]}.parquet", "rows": write_table(table, table_path(directory, name), _table_key(name, table))})
    provenance = dict(provenance or {"package": "sceptre_format.py", "version": "prototype", "language": "Python"})
    provenance.setdefault("written_at", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    _write_json({"stage": stage, **record, "provenance": provenance, "files": files}, _stage_file(directory, stage))


def clear_stage(directory, stage) -> None:
    """Remove a stage and every later stage: their folders, records first (A4).
    The fit cache goes too when the QC stage or an earlier one is cleared, and
    stays when only a testing stage is (T6)."""
    directory = Path(directory)
    cleared = [s for s in STAGES if s == stage or STAGE_RANK[s] > STAGE_RANK[stage]]
    for later in cleared:
        # the record goes first: a stage without a record has not run
        _stage_file(directory, later).unlink(missing_ok=True)
        shutil.rmtree(directory / later, ignore_errors=True)
    if TABLE_SPEC["fits"]["cleared_with"] in cleared:
        shutil.rmtree(directory / FITS, ignore_errors=True)


# The fit cache (T6): fits/response/, a folder of part files and nothing else,
# outside the stages. It is no stage's output and is listed in no record. Each
# writer adds a part of its own, and each part says how its fits were made, in the
# Parquet file's own key-value metadata. No file is shared between writers, so
# there is nothing to keep in step: a part made for other cells, another formula
# or another method is simply not used.
FITS = TABLE_SPEC["fits"]["folder"]
FIT_KEY = TABLE_SPEC["fits"]["identity_key"]


def formula_key(text) -> str:
    """A formula is the same formula however it is spaced (T6.5)."""
    return re.sub(r"\s", "", str(text))


def _identity_value(identity: dict, name):
    value = identity.get(name)
    return formula_key(value) if name == "formula" and value is not None else value


def same_identity(a: dict, b: dict) -> bool:
    """Identities are compared as parsed values, never as text."""
    return all(_identity_value(a, name) == _identity_value(b, name) for name in set(a) | set(b))


def _fit_identity_text(file):
    """The identity a part carries, as the text stored, or None."""
    import pyarrow.parquet as pq

    value = (pq.ParquetFile(file).metadata.metadata or {}).get(FIT_KEY.encode())
    return None if value is None else value.decode()


def fit_parts(directory) -> list:
    """Every part of the cache with the identity it carries: (file, identity). The
    identity is None where a part has none or it is not a JSON object; a reader
    ignores such a part (T6.5)."""
    out = []
    folder = Path(directory) / FITS
    for file in _table_files(folder) if folder.is_dir() else []:
        text = _fit_identity_text(file)
        try:
            identity = None if text is None else json.loads(text)
        except ValueError:
            identity = None
        out.append((file, identity if isinstance(identity, dict) else None))
    return out


def _analysis_fit_fields(directory):
    """What the folder alone says of the identity its fits must have: everything
    but the cells in use, whose hash needs the data file's order of cells."""
    directory = Path(directory)
    if not {"parameters", "qc"} <= set(stages_run(directory)):
        return None
    settings = _read_json(_stage_file(directory, "parameters"))["settings"]
    return {
        "screen_id": _read_json(directory / "analysis.json")["screen"]["id"],
        "formula": settings["formula"],
        "response_fit_method": settings["response_fit_method"],
        "control_group": settings["control_group"],
    }


def _matches(identity, fields) -> bool:
    return identity is not None and all(_identity_value(identity, k) == _identity_value(fields, k) for k in fields)


def write_fits(directory, table: pd.DataFrame, identity: dict) -> Path:
    """Add a part holding `table`, with the identity of its fits inside it."""
    directory = Path(directory)
    if not {"parameters", "qc"} <= set(stages_run(directory)):
        raise ValueError("Fits can be saved once the parameters and qc stages are written.")
    part = directory / FITS / f"part-{uuid.uuid4().hex[:12]}.parquet"
    write_table(table, part, _table_key(FITS, table), {FIT_KEY: json.dumps(identity, sort_keys=True)})
    return part


def read_fits(directory, identity=None):
    """The cached fits made as `identity` says, and that identity, or None if no
    part has it. Parts with another identity, or with none, are ignored. Two parts
    may hold the same response; a reader takes either.

    Without `identity`, the parts made for this analysis's own parameters are
    returned. If those were made for more than one set of cells (parts left by an
    earlier QC), the folder alone cannot say which set is current, and None is
    returned: give the identity."""
    parts = fit_parts(directory)
    if identity is None:
        fields = _analysis_fit_fields(directory)
        if fields is None:
            return None
        found = [i for _, i in parts if _matches(i, fields)]
        if not found or not all(same_identity(i, found[0]) for i in found):
            return None
        identity = found[0]
    files = [file for file, i in parts if i is not None and same_identity(i, identity)]
    if not files:
        return None
    table = read_table(files, TABLES[FITS][2]).drop_duplicates(subset="response_id").reset_index(drop=True)
    return {"identity": identity, "table": table}


def read_analysis(directory) -> dict:
    """An analysis folder: its link, the record of each stage that ran, their
    tables (named as in TABLES), and the fit cache if there is one."""
    directory = Path(directory)
    analysis = _read_json(directory / "analysis.json")
    check_format(analysis.get("format"), directory)
    stages = {s: _read_json(_stage_file(directory, s)) for s in stages_run(directory)}
    tables = {}
    for stage, record in stages.items():
        # a record is the record of its folder's stage, and lists files of that folder only (A2.1, A5.3)
        if record.get("stage") != stage:
            raise ValueError(f"{stage}/record.json says it is the record of the stage '{record.get('stage')}'.")
        for file in record["files"]:
            if not re.fullmatch(r"[A-Za-z0-9_]+\.parquet", file["path"]):
                raise ValueError(f"{stage}/record.json lists '{file['path']}', which is not a table of its own folder.")
            name = f"{stage}/{file['path'].removesuffix('.parquet')}"
            path = table_path(directory, name)
            if not path.exists():
                raise ValueError(f"The stage '{record['stage']}' lists {file['path']}, which is missing.")
            if table_rows(path) != file["rows"]:
                raise ValueError(f"{name}.parquet has {table_rows(path)} rows; the stage '{record['stage']}' wrote {file['rows']}.")
            tables[name] = read_table(path, TABLES[name][2] if name in TABLES else None)
    fits = read_fits(directory)
    if fits is not None:
        tables[FITS] = fits["table"]
    return {
        "format": analysis["format"],
        "response_modality": analysis.get("response_modality"),
        "screen": analysis["screen"],
        "stages": stages,
        "tables": tables,
        "fits_identity": None if fits is None else fits["identity"],
    }


def check_link(recorded: dict, mdata, what) -> None:
    """Refuse a data file other than the one the analysis was run on (A1)."""
    current = screen_link(mdata)
    for field in ("id", "cell_ids_sha256", "response_ids_sha256", "grna_ids_sha256"):
        if recorded.get(field) != current[field]:
            named = f", the dataset '{recorded['name']}'" if recorded.get("name") else ""
            raise ValueError(f"The analysis {what} was run on another data file{named} (its `{field}` differs).")


# ---------------------------------------------------------------- derivations

# Values a package holds that follow from the files; both languages derive them
# the same way.


def top_grna(counts, grna_ids) -> pd.DataFrame:
    """Each cell's gRNA with the most UMIs (the first on a tie) and its share of the
    cell's gRNA UMIs (S6)."""
    from scipy import sparse

    coo = sparse.coo_matrix(counts)
    order = np.lexsort((coo.col, -coo.data.astype(float), coo.row))
    rows, first = np.unique(coo.row[order], return_index=True)
    max_grna = np.full(counts.shape[0], None, dtype=object)
    max_grna[rows] = np.asarray(grna_ids, dtype=object)[coo.col[order][first]]
    fraction = np.full(counts.shape[0], np.nan)
    fraction[rows] = coo.data[order][first] / np.asarray(counts.sum(axis=1)).reshape(-1)[rows]
    return pd.DataFrame({"max_grna": max_grna, "max_grna_frac_umis": fraction})


def grna_groups(grna_table: pd.DataFrame, strategy: str) -> np.ndarray:
    """The gRNA group of each gRNA under an integration strategy."""
    if strategy == "union":
        return grna_table["grna_target"].to_numpy(dtype=object)
    targets = grna_table["grna_target"].to_numpy(dtype=object)
    return np.where(targets == "non-targeting", "non-targeting", grna_table["grna_id"].to_numpy(dtype=object))


def expand_pairs(pairs: pd.DataFrame, grna_table: pd.DataFrame) -> pd.DataFrame:
    """The pairs as a package tests them: one row per (gRNA group, response)."""
    groups = grna_table.loc[grna_table["grna_group"] != "non-targeting", ["grna_target", "grna_group"]].drop_duplicates()
    out = pairs.merge(groups, on="grna_target", how="inner", sort=False)
    return out[["grna_target", "response_id", "grna_group"]].reset_index(drop=True)


def cells_by_grna(assignments: pd.DataFrame, grna_ids, cells: pd.Index) -> dict:
    """Each gRNA's cells (positions among the data file's cells, ascending) (T1)."""
    position = cells.get_indexer(assignments["cell_id"])
    grouped = pd.Series(position).groupby(assignments["grna_id"].to_numpy(), sort=False)
    found = {k: np.sort(v.to_numpy(dtype=np.int64)) for k, v in grouped}
    return {g: found.get(g, np.empty(0, dtype=np.int64)) for g in grna_ids}


def group_cells(initial: dict, grna_table: pd.DataFrame):
    """Each targeting group's cells (gRNA by gRNA, each cell once) and each
    non-targeting gRNA's cells."""
    targeting = (grna_table["grna_group"] != "non-targeting").to_numpy()
    groups = {}
    for group in dict.fromkeys(grna_table.loc[targeting, "grna_group"]):
        ids = grna_table.loc[grna_table["grna_group"] == group, "grna_id"]
        groups[group] = pd.unique(np.concatenate([initial[i] for i in ids])).astype(np.int64)
    return groups, {i: initial[i] for i in grna_table.loc[~targeting, "grna_id"]}


def grnas_per_cell(initial: dict, n_cells: int) -> np.ndarray:
    return np.bincount(np.concatenate(list(initial.values())).astype(np.int64), minlength=n_cells)


def not_one_grna(initial: dict, n_cells: int) -> np.ndarray:
    """In low MOI, the cells that do not carry exactly one gRNA."""
    return np.flatnonzero(grnas_per_cell(initial, n_cells) != 1)


def qc_assignments(groups: dict, nt: dict, cells_in_use: np.ndarray, nt_cells_control: bool):
    """After QC: positions among the cells in use; with the NT cells as the control
    group, each non-targeting gRNA's cells as positions among all NT cells."""
    size = int(cells_in_use.max()) + 1 if cells_in_use.size else 0
    position = np.full(size, -1, dtype=np.int64)
    position[cells_in_use] = np.arange(cells_in_use.size)

    def among(cells):
        found = position[cells[cells < size]]
        return np.sort(found[found >= 0])

    groups = {k: among(v) for k, v in groups.items()}
    nt = {k: among(v) for k, v in nt.items()}
    all_nt = None
    if nt_cells_control:
        all_nt = pd.unique(np.concatenate(list(nt.values()))).astype(np.int64) if nt else np.empty(0, dtype=np.int64)
        stop = np.cumsum([len(v) for v in nt.values()])
        start = np.concatenate([[0], stop[:-1]])
        nt = {k: np.arange(a, b, dtype=np.int64) for k, a, b in zip(nt, start, stop, strict=True)}
    return groups, nt, all_nt


def n_pairs_passing(pair_qc: pd.DataFrame, strategy: str) -> int:
    """How many pairs pass pairwise QC, as sceptre counts them: under the bonferroni
    strategy a (target, response) pair passes if any of its gRNAs does."""
    passed = pair_qc["pass_qc"].astype(bool)
    if strategy != "bonferroni":
        return int(passed.sum())
    return int(passed.groupby([pair_qc["grna_target"], pair_qc["response_id"]]).any().sum())


def significance_threshold(result: pd.DataFrame):
    """The largest p-value called significant, or None: `significant` is then
    p_value <= threshold, which no part file needs the whole table for (T5)."""
    called = result["significant"].fillna(False).to_numpy(dtype=bool)
    return float(result.loc[called, "p_value"].max()) if called.any() else None


def add_significant(result: pd.DataFrame, threshold) -> pd.DataFrame:
    result = result.copy()
    missing = result["p_value"].isna()
    called = pd.array(np.zeros(len(result), dtype=bool) if threshold is None else (result["p_value"] <= threshold).to_numpy(), dtype="boolean")
    called[missing.to_numpy()] = pd.NA
    result["significant"] = called
    return result


def log_2_of(fold_change) -> np.ndarray:
    """sceptre's formula, to the last bit: log(x) / log(2), which is not always what
    log2(x) returns."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(np.asarray(fold_change, dtype=float)) / np.log(2.0)


def log_2_follows(result: pd.DataFrame) -> bool:
    """Whether `log_2_fold_change` follows from `fold_change`; a results table then
    does not store it (T5)."""
    if "fold_change" not in result or "log_2_fold_change" not in result:
        return False
    derived = log_2_of(result["fold_change"])
    held = result["log_2_fold_change"].to_numpy(dtype=float)
    return bool(np.all((np.isnan(derived) & np.isnan(held)) | (derived == held)))


def add_log_2_fold_change(result: pd.DataFrame) -> pd.DataFrame:
    if "log_2_fold_change" in result or "fold_change" not in result:
        return result
    result = result.copy()
    after = max(list(result.columns).index(c) for c in ("fold_change", "se_fold_change") if c in result)
    result.insert(after + 1, "log_2_fold_change", log_2_of(result["fold_change"]))
    return result


def resample_columns(table: pd.DataFrame) -> list[str]:
    return [c for c in table.columns if re.fullmatch(r"z_null_[0-9]+", str(c))]


def resamples_table(result: pd.DataFrame):
    """Resampled statistics: packages hold them as columns z_null_1, z_null_2, ...
    padded with missing values; the file holds one list per pair (T5)."""
    columns = resample_columns(result)
    if not columns:
        return None
    values = result[columns].to_numpy(dtype=float)
    out = result[unit_key(result)].reset_index(drop=True).copy()
    out["z_null"] = [row[~np.isnan(row)].tolist() for row in values]
    return out


def add_resamples(result: pd.DataFrame, resamples) -> pd.DataFrame:
    if resamples is None:
        return result
    key = unit_key(result)
    lists = result[key].merge(resamples, on=key, how="left")["z_null"]
    width = max(len(z) for z in lists)
    values = np.full((len(result), width), np.nan)
    for i, z in enumerate(lists):
        values[i, : len(z)] = z
    wide = pd.DataFrame(values, columns=[f"z_null_{j + 1}" for j in range(width)], index=result.index)
    return pd.concat([result, wide], axis=1)


def formula_in_subset(text: str) -> bool:
    """Additive formulas only: each term a covariate name, log(x) or log(x + 1),
    with 0 or 1 for the intercept (T7)."""
    rhs = re.sub(r"\s", "", text.split("~", 1)[1])
    terms, depth, start = [], 0, 0
    for i, char in enumerate(rhs):
        depth += (char == "(") - (char == ")")
        if char == "+" and depth == 0:
            terms.append(rhs[start:i])
            start = i + 1
    terms.append(rhs[start:])
    name = "[A-Za-z][A-Za-z0-9_.]*"
    return all(re.fullmatch(rf"0|1|{name}|log\({name}\)|log\({name}\+1\)", t) for t in terms)


def design_sha256(matrix) -> str:
    """SHA-256 of the matrix's values as little-endian doubles, column by column."""
    return hashlib.sha256(np.asarray(matrix, dtype="<f8").tobytes(order="F")).hexdigest()


# ---------------------------------------------------------------- comparing and describing


def rows_sorted(table: pd.DataFrame) -> pd.DataFrame:
    """A table with its columns and rows in a fixed order, for comparing."""
    columns = sorted(map(str, table.columns))
    plain = table.rename(columns=str)[columns].astype(object)
    plain = plain.where(plain.notna(), None)
    order = sorted(range(len(plain)), key=lambda i: tuple(repr(v) for v in plain.iloc[i]))
    return plain.iloc[order].reset_index(drop=True)


def same_table(a: pd.DataFrame, b: pd.DataFrame) -> bool:
    """Whether two tables hold the same rows, in any order."""
    if len(a) != len(b) or sorted(map(str, a.columns)) != sorted(map(str, b.columns)):
        return False
    x, y = rows_sorted(a), rows_sorted(b)
    for column in x.columns:
        for u, v in zip(x[column], y[column], strict=True):
            if isinstance(u, (list, np.ndarray)) or isinstance(v, (list, np.ndarray)):
                if not np.array_equal(np.asarray(u, dtype=float), np.asarray(v, dtype=float)):
                    return False
            elif u != v and not (u is None and v is None):
                return False
    return True


def file_listing(directory) -> pd.DataFrame:
    directory = Path(directory)
    files = sorted(p for p in directory.rglob("*") if p.is_file())
    return pd.DataFrame({"file": [str(p.relative_to(directory)) for p in files], "kilobytes": [round(p.stat().st_size / 1024, 1) for p in files]})


# ---------------------------------------------------------------- validation

# validate_screen() and validate_analysis() check files against the numbered rules
# of SPEC.md. They read the files directly, build no object, and report every
# problem they find: a DataFrame with one row per problem (`rule`, `where`,
# `message`), with no rows when the files conform. The R library has the same two
# functions and reports the same rules.

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
            log.add("S10.1", where, "needs the scalar text entries `format`, `moi` and `id`.")
        elif format_parts(_text(record["format"][()])[0]) is None:
            log.add("S10.1", f"{where}/format", "is not a version of the form MAJOR.MINOR.")
        else:
            if _text(record["moi"][()])[0] not in entries["moi"]["values"]:
                log.add("S10.2", f"{where}/moi", "must be 'high' or 'low'.")
            if not known_major(_text(record["format"][()])[0]):
                log.add("S10.3", f"{where}/format", "is a format this validator does not know.")
            if not _text(record["id"][()])[0]:
                log.add("S10.1", f"{where}/id", "is empty.")
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
                    log.add("S5.4", path, "is not a valid covariate name (letters, digits, underscores and dots, starting with a letter).")
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
    return log.all()


def _screen_identity(file, response_modality) -> dict:
    """What validate_analysis() needs from a data file, read without building an object."""
    import h5py

    with h5py.File(file, "r") as f:
        out = {
            "modalities": list(f["mod"]),
            "id": _text(f["uns/sceptre/id"][()])[0],
            "cells": _index_of(f["mod/grna/obs"]),
            "grnas": _index_of(f["mod/grna/var"]),
        }
        target = f["mod/grna/var/grna_target"]
        if _encoding(target) == "categorical":
            categories = _text(target["categories"][()])
            out["targets"] = [categories[c] for c in target["codes"][()]]
        else:
            out["targets"] = _text(target[()])
        if response_modality in f["mod"]:
            out["responses"] = _index_of(f[f"mod/{response_modality}/var"])
    return out


def _schema_problems(path, schema) -> list[str]:
    """Problems a JSON file has against a schema, as text; empty when it conforms."""
    return _schema_problems_text(Path(path).read_text(), schema)


def _schema_problems_text(text, schema) -> list[str]:
    import jsonschema

    try:
        value = json.loads(text)
    except ValueError:
        return ["is not valid JSON."]
    validator = jsonschema.Draft7Validator(json.loads((SPEC_DIR / schema).read_text()))
    return sorted({("/" + "/".join(map(str, e.absolute_path))) + " " + e.message for e in validator.iter_errors(value)})


def _column_kind(type_) -> str:
    import pyarrow as pa

    if pa.types.is_dictionary(type_):
        type_ = type_.value_type
    if pa.types.is_string(type_) or pa.types.is_large_string(type_):
        return "string"
    if pa.types.is_boolean(type_):
        return "boolean"
    if pa.types.is_floating(type_):
        return "double"
    if (pa.types.is_list(type_) or pa.types.is_large_list(type_)) and pa.types.is_floating(type_.value_type):
        return "list<double>"
    return "other"


def validate_analysis(directory, screen_file=None) -> pd.DataFrame:
    import pyarrow.parquet as pq

    directory = Path(directory)
    log = _Problems()
    path = directory / "analysis.json"
    if not path.exists():
        log.add("A1.1", "analysis.json", "is missing.")
        return log.all()
    faults = _schema_problems(path, "analysis.schema.json")
    for fault in faults:
        # the form of the hint has a rule of its own
        log.add("A1.4" if fault.startswith("/screen/path ") else "A1.1", "analysis.json", fault)
    analysis = None if faults else _read_json(path)
    if analysis is not None and not known_major(analysis["format"]):
        log.add("A1.1", "analysis.json", "is in a format this validator does not know.")

    # stage records: <stage>/record.json
    present = stages_run(directory)
    records, unchecked = {}, []   # unchecked: stages whose record cannot be trusted
    for stage in present:
        where = f"{stage}/record.json"
        faults = _schema_problems(_stage_file(directory, stage), "stage.schema.json")
        for fault in faults:
            log.add("A2.1", where, fault)
        if faults:
            unchecked.append(stage)
            continue
        record = _read_json(_stage_file(directory, stage))
        if record["stage"] != stage:
            log.add("A2.1", where, f"says it is the record of the stage '{record['stage']}'.")
            unchecked.append(stage)
            continue
        records[stage] = record
    for stage in present:
        missing = [s for s in STAGES if STAGE_RANK[s] < STAGE_RANK[stage] and s not in present]
        if missing:
            log.add("A2.2", f"{stage}/record.json", f"The stage needs the stage(s) {', '.join(missing)} first.")

    # the files each record lists, all in the stage's own folder
    listed, tables = [], {}
    for stage, record in records.items():
        where = f"{stage}/record.json"
        names_listed = [f"{stage}/{f['path'].removesuffix('.parquet')}" for f in record["files"]]
        listed += names_listed
        for file, name in zip(record["files"], names_listed, strict=True):
            shown, full = f"{name}.parquet", table_path(directory, name)
            if not full.exists():
                log.add("A5.1", shown, "is listed by its stage's record and is missing.")
                continue
            try:
                rows = table_rows(full)
            except Exception:
                log.add("A5.1", shown, "cannot be read as Parquet.")
                continue
            if rows != file["rows"]:
                log.add("A5.2", shown, f"has {rows} rows; the stage '{stage}' wrote {file['rows']}.")
            if name in TABLE_DEFINITIONS:
                tables[name] = full
        for name, definition in TABLE_DEFINITIONS.items():
            if definition["stage"] == stage and definition["required"] and name not in names_listed:
                log.add("A5.4", where, f"does not list the table {name}.")
    # a stage's table that is in its folder and that no record lists (A4.1)
    for name, definition in TABLE_DEFINITIONS.items():
        stage = definition["stage"]
        if stage is None or stage in unchecked or name in listed:
            continue
        if table_path(directory, name).exists():
            log.add("A4.1", f"{name}.parquet", "is in the folder and its stage's record does not list it.")

    # the fit cache: part files only, each saying how its fits were made (T6)
    fit_identities = {}
    if (directory / FITS).is_dir():
        for entry in sorted((directory / FITS).iterdir()):
            if entry.is_dir() or entry.suffix != ".parquet":
                log.add("T6.1", f"{FITS}/{entry.name}", "is not a part file; the fit cache holds part files only.")
        for file in _table_files(directory / FITS):
            try:
                text = _fit_identity_text(file)
            except Exception:
                continue   # not a Parquet file: reported with the table below
            if text is None:
                log.add("T6.1", f"{FITS}/{file.name}", f"does not say how its fits were made: its metadata has no `{FIT_KEY}`.")
                continue
            faults = _schema_problems_text(text, TABLE_SPEC["fits"]["identity_schema"])
            for fault in faults:
                log.add("T6.1", f"{FITS}/{file.name}", f"`{FIT_KEY}` {fault}")
            if not faults:
                fit_identities[file] = json.loads(text)
        if _table_files(directory / FITS):
            tables[FITS] = directory / FITS
        if not {"parameters", "qc"} <= set(present):
            log.add("T6.3", FITS, "The fit cache needs the parameters and qc stages.")

    # each table against its definition
    strategy = records.get("parameters", {}).get("settings", {}).get("grna_integration_strategy")
    contents = {}
    for name, full in tables.items():
        definition, where = TABLE_DEFINITIONS[name], (name if name == FITS else f"{name}.parquet")
        parts = []
        for file in _table_files(full):
            parquet = pq.ParquetFile(file)
            if any(parquet.metadata.row_group(i).num_rows > ROW_GROUP_ROWS for i in range(parquet.metadata.num_row_groups)):
                log.add("T0.3", where, f"has a row group of more than {ROW_GROUP_ROWS} rows.")
            parts.append(parquet.read())
        ok = True
        for column in definition["columns"]:
            for part in parts:
                if column["name"] not in part.column_names:
                    if column["required"]:
                        log.add("T0.1", where, f"has no column `{column['name']}`.")
                        ok = False
                elif _column_kind(part.schema.field(column["name"]).type) != column["type"]:
                    log.add("T0.1", where, f"The column `{column['name']}` is not of type {column['type']}.")
                    ok = False
                elif not column.get("nullable") and column["type"] != "list<double>" and part[column["name"]].null_count:
                    log.add("T0.1", where, f"The column `{column['name']}` has missing values.")
                    ok = False
        for group in definition.get("one_of", []):
            if not any(c in parts[0].column_names for c in group):
                log.add("T5.1", where, f"needs one of the columns {', '.join(group)}.")
        for column in definition.get("forbidden_columns", []):
            if column["name"] in parts[0].column_names:
                log.add(column["rule"], where, f"must not have the column `{column['name']}`.")
        # a column present under some gRNA integration strategies only (T3.3, T5.6)
        for column in definition["columns"]:
            if "present_under" not in column or strategy is None:
                continue
            expected = strategy in column["present_under"]
            if expected != (column["name"] in parts[0].column_names):
                verb = "needs" if expected else "must not have"
                log.add(column["rule"], where, f"{verb} the column `{column['name']}` under the strategy '{strategy}'.")
        if not ok:
            continue
        # within one file the key is sorted and unique; across the parts of the fit
        # cache it need not be unique
        across_parts = definition.get("unique_across_parts", True)
        frames = [read_table_from_arrow(part) for part in parts]
        for frame in frames:
            key = _table_key(name, frame)
            if len(frame) > 1 and _byte_order(frame, key) != list(range(len(frame))):
                log.add("T0.2", where, f"The rows are not sorted by {', '.join(key)}.")
            if not across_parts and frame.duplicated(subset=key).any():
                log.add("T0.4", where, f"Two rows of one file have the same {', '.join(key)}.")
        # the parts already read, put together: each table is read once
        whole = frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True, join="inner")
        key = _table_key(name, whole)
        if across_parts and whole.duplicated(subset=key).any():
            log.add("T0.4", where, f"Two rows have the same {', '.join(key)}.")
        contents[name] = whole

    parameters = records.get("parameters")
    if parameters is not None:
        if not formula_in_subset(parameters["settings"]["formula"]):
            log.add("T7.1", "parameters/record.json", "The formula is outside the additive subset.")
        # A part made for this analysis's parameters has one coefficient per term of its
        # design (T6.2). A part made otherwise is allowed and unused, and is not
        # reported: nothing depends on it, and clearing it is housekeeping (T6.3, T6.6).
        expected = {
            "formula": parameters["settings"]["formula"],
            "response_fit_method": parameters["settings"]["response_fit_method"],
            "control_group": parameters["settings"]["control_group"],
        }
        if analysis is not None:
            expected["screen_id"] = analysis["screen"]["id"]
        for file, identity in fit_identities.items():
            if _matches(identity, expected) and "coefficients" in pq.ParquetFile(file).schema_arrow.names:
                if any(len(c) != len(parameters["design"]["terms"]) for c in pq.read_table(file, columns=["coefficients"])["coefficients"].to_pylist() if c is not None):
                    log.add("T6.2", FITS, "A fit does not have one coefficient per term of the design.")
                    break

    # against the data file
    if screen_file is not None and analysis is not None:
        screen = _screen_identity(screen_file, analysis["response_modality"])
        if "responses" not in screen:
            log.add("A1.2", "analysis.json", f"The data file has no modality named '{analysis['response_modality']}'.")
        else:
            current = {
                "id": screen["id"], "n_cells": len(screen["cells"]), "n_responses": len(screen["responses"]), "n_grnas": len(screen["grnas"]),
                "cell_ids_sha256": sha256_ids(screen["cells"]), "response_ids_sha256": sha256_ids(screen["responses"]),
                "grna_ids_sha256": sha256_ids(screen["grnas"]),
            }
            for field, value in current.items():
                if analysis["screen"].get(field) != value:
                    log.add("A1.2", "analysis.json", f"`screen.{field}` is not that of the data file.")
            known = {"cell_id": set(screen["cells"]), "response_id": set(screen["responses"])}
            for name, table in contents.items():
                for column, values in known.items():
                    if column in table and not set(table[column]) <= values:
                        log.add("T0.5", f"{name}.parquet", f"`{column}` names something the data file does not have.")
            for name in ("assignment/assignments", "calibration/negative_control_groups"):
                if name in contents and not set(contents[name]["grna_id"]) <= set(screen["grnas"]):
                    log.add("T0.5", f"{name}.parquet", "`grna_id` names something the data file does not have.")
            for name in ("parameters/discovery_pairs", "parameters/positive_control_pairs"):
                if name in contents and not set(contents[name]["grna_target"]) <= set(screen["targets"]):
                    log.add("T0.5", f"{name}.parquet", "`grna_target` names something the data file does not have.")
            kept = contents.get("qc/cells")
            if kept is not None and set(kept["cell_id"]) != set(screen["cells"]):
                log.add("T2.1", "qc/cells.parquet", "does not have one row per cell of the data file.")
    return log.all()
